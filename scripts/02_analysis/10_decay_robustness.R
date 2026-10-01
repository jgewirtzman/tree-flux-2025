# ============================================================
# 10_decay_robustness.R
#
# How robust are the specialist decay-flux correlations to (i) how the fluxes
# were calculated and (ii) how each tree's flux is summarized?
#
# Flux versions (tree means over the same closures):
#   legacy        fluxes in the published table / bioRxiv v1 (upstream linear fits)
#   goflux_dbXX   goFlux refits with a 0/10/20/30/45-s deadband (13_deadband_sensitivity.R)
#   current       canonical fluxes of the analysis dataset (20-s deadband, end trim, manual review)
# Tree summaries: mean of all measurements, median, mean of detected fluxes only (|flux| > MDF),
#   mean of growing-season (May-Oct) measurements.
# Metrics: ERT CV, ERT mean resistivity, ERT index (species-normalized PC1, tomography paper),
#   SoT structural loss.
#
# Outputs: outputs/tables/decay_robustness_correlations.csv  (species x version x summary x metric)
#          outputs/tables/decay_robustness_tree_means.csv    (per-tree values, Q. rubra and N. sylvatica)
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
g <- read.csv("data/input/HF_2023-2025_tree_flux_goflux.csv", stringsAsFactors = FALSE)
an <- read.csv("data/processed/flux_with_quality_flags.csv", stringsAsFactors = FALSE)
db <- read.csv("outputs/tables/deadband_sensitivity.csv", stringsAsFactors = FALSE)
tc <- read.csv("data/processed/tomography_classes.csv", stringsAsFactors = FALSE)
ert <- read.csv("data/input/tomography_results_compiled.csv", stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
met <- ert %>% transmute(Tree = tree, ert_cv, ert_mean) %>%
  inner_join(tc %>% transmute(Tree = tree, ert_pc1, sot_loss = sot_structural_loss), by = "Tree")

keep_ids <- an$closure_id                                # closures in the analysis dataset
base <- an %>% select(closure_id, current = CH4_flux_nmolpm2ps, CH4_MDF) %>%
  inner_join(g %>% select(closure_id, Tree, SPECIES, location, date, legacy = CH4_flux_legacy), by = "closure_id")
dbw <- db %>% filter(closure_id %in% keep_ids) %>% select(closure_id, deadband_s, flux) %>%
  mutate(v = sprintf("goflux_db%02d", deadband_s)) %>% select(-deadband_s) %>% pivot_wider(names_from = v, values_from = flux)
d <- base %>% left_join(dbw, by = "closure_id") %>%
  mutate(month = as.integer(substr(date, 6, 7)))
versions <- c("legacy", grep("^goflux_db", names(d), value = TRUE), "current")
long <- d %>% pivot_longer(all_of(versions), names_to = "version", values_to = "flux") %>% filter(!is.na(flux))
summ <- bind_rows(
  long %>% group_by(Tree, SPECIES, location, version) %>% summarise(value = mean(flux), .groups = "drop") %>% mutate(summary = "mean"),
  long %>% group_by(Tree, SPECIES, location, version) %>% summarise(value = median(flux), .groups = "drop") %>% mutate(summary = "median"),
  long %>% filter(abs(flux) > CH4_MDF) %>% group_by(Tree, SPECIES, location, version) %>% summarise(value = mean(flux), .groups = "drop") %>% mutate(summary = "mean_detected"),
  long %>% filter(month %in% 5:10) %>% group_by(Tree, SPECIES, location, version) %>% summarise(value = mean(flux), .groups = "drop") %>% mutate(summary = "mean_growing_season"))
x <- summ %>% inner_join(met, by = "Tree")
res <- x %>% pivot_longer(c(ert_cv, ert_mean, ert_pc1, sot_loss), names_to = "metric", values_to = "m") %>%
  filter(!is.na(m)) %>% group_by(SPECIES, location, version, summary, metric) %>%
  summarise(n = n(), r = cor(m, value), p = cor.test(m, value)$p.value,
            rho = suppressWarnings(cor(m, value, method = "spearman")),
            rho_p = suppressWarnings(cor.test(m, value, method = "spearman")$p.value),
            loo_min = min(sapply(seq_len(n()), function(i) cor(m[-i], value[-i]))),
            loo_max = max(sapply(seq_len(n()), function(i) cor(m[-i], value[-i]))), .groups = "drop")
write.csv(res, "outputs/tables/decay_robustness_correlations.csv", row.names = FALSE)
tm <- x %>% filter(SPECIES %in% c("ro", "bg"), summary == "mean") %>% select(Tree, SPECIES, version, value, ert_cv, ert_pc1, ert_mean, sot_loss) %>%
  pivot_wider(names_from = version, values_from = value)
write.csv(tm, "outputs/tables/decay_robustness_tree_means.csv", row.names = FALSE)
options(width = 220)
show <- function(sp) {
  cat("\n=====", sp, "— Pearson r (p) of tree-level flux with each metric =====\n")
  y <- res %>% filter(SPECIES == sp) %>% mutate(cell = sprintf("%5.2f (%.3f)", r, p)) %>%
    select(summary, version, metric, cell) %>% pivot_wider(names_from = metric, values_from = cell) %>% arrange(summary, version)
  print(as.data.frame(y), row.names = FALSE)
}
show("ro"); show("bg")
cat("\n===== Q. rubra tree means by flux version =====\n")
print(as.data.frame(tm %>% filter(SPECIES == "ro") %>% arrange(ert_pc1) %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)

# ------------------------------------------------------------
# Figure: tree-mean flux (current) vs each wood-condition metric, both specialists.
# Line drawn only where Pearson p < 0.05; subtitle gives r, rho and r without the
# highest-index tree (the single most influential point for Q. rubra).
# ------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
lab <- c(ert_cv = "ERT CV", ert_mean = "ERT mean (Ohm m)", ert_pc1 = "ERT index (PC1)", sot_loss = "SoT structural loss (%)")
spl <- c(bg = "N. sylvatica (wetland)", ro = "Q. rubra (upland)")
pd <- x %>% filter(SPECIES %in% names(spl), summary == "mean", version == "current") %>%
  pivot_longer(names(lab), names_to = "metric", values_to = "m") %>% filter(!is.na(m)) %>%
  mutate(metric = factor(lab[metric], lab), sp = factor(spl[SPECIES], spl))
st <- pd %>% group_by(sp, metric) %>% summarise(
  r = cor(m, value), p = cor.test(m, value)$p.value,
  rho = suppressWarnings(cor(m, value, method = "spearman")),
  r_drop = { i <- which.max(m); cor(m[-i], value[-i]) }, .groups = "drop") %>%
  mutate(txt = sprintf("r = %.2f (p = %.3f)%s\nrho = %.2f; r without top tree = %.2f", r, p, ifelse(p < 0.05, " *", ""), rho, r_drop))
pfig <- ggplot(pd, aes(m, value)) +
  geom_smooth(data = pd %>% semi_join(st %>% filter(p < 0.05), by = c("sp", "metric")),
              method = "lm", formula = y ~ x, colour = "#2a7f74", fill = "#2a7f74", alpha = .12) +
  geom_point(size = 2.4, shape = 21, fill = "grey35", colour = "white") +
  geom_text(data = st, aes(x = -Inf, y = Inf, label = txt), hjust = -0.03, vjust = 1.15, size = 2.9, inherit.aes = FALSE) +
  facet_grid(sp ~ metric, scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(.05, .35))) +
  labs(x = NULL, y = expression("Tree-mean CH"[4]*" flux (nmol m"^-2*" s"^-1*")")) +
  theme_bw(base_size = 11) + theme(strip.background = element_rect(fill = "grey95"), panel.grid.minor = element_blank())
ggsave("outputs/figures/tomography/decay_metric_correlations_specialists.png", pfig, width = 12, height = 6.2, dpi = 200)
