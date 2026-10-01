# ============================================================
# 07_decay_definitions.R
#
# Does the decay-flux result depend on how internal wood condition is defined?
# Compares, for each site (and species within site):
#   ert_cv          raw ERT coefficient of variation (bioRxiv v1)
#   ert_cv_z        ERT CV z-scored within species
#   ert_pc1         species-normalized ERT PC1 (tomography paper; classification axis)
#   moisture_anom   PC1 > study-set mean (binary, tomography-paper threshold)
#   sot_loss        SoT structural loss (%), sqrt-transformed
#   decay_any       class II-IV vs I (tomography-paper classes)
#   ert_mean_z      mean resistivity z-scored within species (low = wet)
#   ert_cma         central moisture accumulation
# Three analyses per definition:
#   (1) Pearson r and Spearman rho with tree-mean CH4 flux (n = 30 trees per site)
#   (2) leave-one-tree-out range of r (influence of single trees)
#   (3) species-adjusted mixed model on ALL observations:
#       asinh(1000 x flux) ~ metric + species + (1 | tree), per site
#       (species-adjusted = within-species association; avoids confounding species
#        differences in flux with species differences in wood properties)
#
# Output: outputs/tables/decay_definition_comparison.csv
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(lme4) })

tc <- read.csv("data/final/tomography_classes.csv", stringsAsFactors = FALSE)
ert <- read.csv("data/package/tomography/tomography_results_compiled.csv", stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
fx <- read.csv("data/final/stem_ch4_flux.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(PLOT), !is.na(CH4_flux_nmolpm2ps)) %>%
  group_by(Tree) %>% mutate(SPECIES = first(na.omit(SPECIES))) %>% ungroup() %>%
  mutate(site = ifelse(PLOT == "BGS", "Wetland", "Upland"), y = asinh(1000 * CH4_flux_nmolpm2ps))

zs <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
trees <- ert %>% transmute(tree = as.character(tree), ert_cv, ert_mean, ert_cma) %>%
  inner_join(tc %>% mutate(tree = as.character(tree)), by = "tree") %>%
  inner_join(fx %>% group_by(tree = as.character(Tree)) %>%
               summarise(flux = mean(CH4_flux_nmolpm2ps), SPECIES = first(SPECIES), .groups = "drop"), by = "tree") %>%
  mutate(ert_cv_z = zs(ert_cv, SPECIES), ert_mean_z = zs(ert_mean, SPECIES), ert_pc1 = ert_pc1,
         moisture_anom = as.numeric(ert_pc1 > ert_pc1_threshold),
         sot_loss = sqrt(sot_structural_loss), decay_any = as.numeric(decay_class != "I: No Decay"))

defs <- c("ert_cv", "ert_cv_z", "ert_pc1", "moisture_anom", "sot_loss", "decay_any", "ert_mean_z", "ert_cma")
res <- list()
for (s in c("Wetland", "Upland")) for (v in defs) {
  d <- trees %>% filter(site == s, !is.na(.data[[v]]))
  ct <- cor.test(d[[v]], d$flux); sp <- suppressWarnings(cor.test(d[[v]], d$flux, method = "spearman"))
  loo <- vapply(seq_len(nrow(d)), function(i) cor(d[[v]][-i], d$flux[-i]), numeric(1))
  o <- fx %>% filter(site == s) %>% mutate(tree = as.character(Tree)) %>%
    inner_join(d %>% select(tree, m = all_of(v)), by = "tree")
  o$m <- as.numeric(scale(o$m))
  mm <- lmer(y ~ m + SPECIES + (1 | Tree), data = o, REML = FALSE)
  m0 <- update(mm, . ~ . - m)
  co <- summary(mm)$coefficients["m", ]
  res[[length(res) + 1]] <- tibble(site = s, definition = v, n_trees = nrow(d),
    pearson_r = ct$estimate, pearson_p = ct$p.value, spearman_rho = sp$estimate, spearman_p = sp$p.value,
    loo_r_min = min(loo), loo_r_max = max(loo),
    mixed_beta_per_sd = co[["Estimate"]], mixed_se = co[["Std. Error"]],
    mixed_lrt_p = anova(m0, mm)$`Pr(>Chisq)`[2], n_obs = nrow(o))
}
res <- bind_rows(res)
# within-species correlations for the specialists (tree means)
sp_res <- list()
for (spp in c("bg", "ro", "rm", "hem")) for (s in c("Wetland", "Upland")) for (v in c("ert_cv", "ert_cv_z", "ert_pc1", "sot_loss", "ert_mean_z")) {
  d <- trees %>% filter(SPECIES == spp, site == s, !is.na(.data[[v]]))
  if (nrow(d) < 6) next
  ct <- cor.test(d[[v]], d$flux); spc <- suppressWarnings(cor.test(d[[v]], d$flux, method = "spearman"))
  sp_res[[length(sp_res) + 1]] <- tibble(site = s, species = spp, definition = v, n = nrow(d),
    r = ct$estimate, p = ct$p.value, rho = spc$estimate, rho_p = spc$p.value)
}
sp_res <- bind_rows(sp_res)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "outputs/tables/decay_definition_comparison.csv", row.names = FALSE)
write.csv(sp_res, "outputs/tables/decay_definition_by_species.csv", row.names = FALSE)
options(width = 200)
print(as.data.frame(res %>% mutate(across(where(is.numeric), ~ signif(.x, 2)))))
print(as.data.frame(sp_res %>% mutate(across(where(is.numeric), ~ signif(.x, 2)))))

# ------------------------------------------------------------
# SI table: correlation of tree-mean CH4 flux with four wood-condition metrics
#   SoT structural loss  % of the SoT cross-section in non-brown (low-velocity) classes (PiCUS Q74)
#   ERT mean             mean resistivity of the ERT cross-section (Ohm m; lower = wetter)
#   ERT CV               coefficient of variation of resistivity (heterogeneity of moisture)
#   ERT index (PC1)      first principal component of eight ERT metrics, each z-scored within
#                        species (tomography paper; higher = wetter, more heterogeneous)
# rows: each species x site (10 trees), each site pooled (30 trees), and the species-adjusted
# mixed-model test on all measurements (p only)
# ------------------------------------------------------------
metrics <- c(sot_loss_pct = "SoT structural loss (%)", ert_mean = "ERT mean (Ohm m)", ert_cv = "ERT CV", ert_pc1 = "ERT index (PC1)")
trees$sot_loss_pct <- trees$sot_structural_loss
fmt <- function(r, p) ifelse(is.na(r), "–", sprintf("%.2f (%s)%s", r, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
                                                      ifelse(p < 0.05, "*", "")))
# Pearson r (p), Spearman rho, and the range of r when each tree is omitted in turn (robustness)
cell <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  ct <- cor.test(x, y); rho <- suppressWarnings(cor(x, y, method = "spearman"))
  loo <- vapply(seq_along(x), function(i) cor(x[-i], y[-i]), numeric(1))
  sprintf("%s; ρ = %.2f; LOO %.2f to %.2f", fmt(ct$estimate, ct$p.value), rho, min(loo), max(loo))
}
grp <- list()
for (s in c("Wetland", "Upland")) {
  for (spp in c("bg", "rm", "hem", "ro")) {
    d <- trees %>% filter(site == s, SPECIES == spp); if (nrow(d) < 6) next
    grp[[length(grp) + 1]] <- c(group = sprintf("%s — %s", s, c(bg = "N. sylvatica", rm = "A. rubrum", hem = "T. canadensis", ro = "Q. rubra")[[spp]]),
      n = nrow(d), sapply(names(metrics), function(m) cell(d[[m]], d$flux)))
  }
  d <- trees %>% filter(site == s)
  grp[[length(grp) + 1]] <- c(group = sprintf("%s — all trees pooled", s), n = nrow(d),
    sapply(names(metrics), function(m) cell(d[[m]], d$flux)))
  # species-adjusted mixed model on all measurements (p of the metric term)
  pv <- sapply(names(metrics), function(m) {
    o <- fx %>% filter(site == s) %>% mutate(tree = as.character(Tree)) %>% inner_join(d %>% select(tree, x = all_of(m)), by = "tree") %>% filter(!is.na(x))
    o$x <- as.numeric(scale(o$x)); mm <- lmer(y ~ x + SPECIES + (1 | Tree), data = o, REML = FALSE)
    b <- fixef(mm)[["x"]]; p <- anova(update(mm, . ~ . - x), mm)$`Pr(>Chisq)`[2]
    sprintf("β = %.2f (%s)%s", b, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)), ifelse(p < 0.05, "*", "")) })
  grp[[length(grp) + 1]] <- c(group = sprintf("%s — species-adjusted, all measurements", s), n = sum(fx$site == s), pv)
}
si <- as.data.frame(do.call(rbind, grp), stringsAsFactors = FALSE)
names(si) <- c("Group", "n", unname(metrics))
write.csv(si, "outputs/tables/SI_decay_metric_correlations.csv", row.names = FALSE)
print(si)
