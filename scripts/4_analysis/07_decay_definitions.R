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
#       asinh(flux) ~ metric + species + (1 | tree), per site (the scale of all other models)
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
  mutate(site = ifelse(PLOT == "BGS", "Wetland", "Upland"), y = asinh(CH4_flux_nmolpm2ps))

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
# SI table: four wood-condition metrics against every CH4 measurement
#   SoT structural loss  % of the SoT cross-section in non-brown (low-velocity) classes (PiCUS Q74)
#   ERT mean             mean resistivity of the ERT cross-section (Ohm m; lower = wetter)
#   ERT CV               coefficient of variation of resistivity (heterogeneity of moisture)
#   ERT index (PC1)      first principal component of eight ERT metrics, each z-scored within
#                        species (tomography paper; higher = wetter, more heterogeneous)
# Each metric is measured once per tree, so the tree is the unit of replication: every
# measurement enters asinh(flux) ~ metric (per SD) + (1 | tree) + (1 | date), tested with
# Kenward-Roger degrees of freedom (about 8 for 10 trees), as in Figure 5 (06_tomography_flux.R).
# Rows: each species x site (10 trees), and per site all species with species as a fixed effect.
# Cells: beta per SD (p), and the range of beta when each tree is omitted in turn.
# ------------------------------------------------------------
suppressPackageStartupMessages(library(lmerTest))
metrics <- c(sot_loss_pct = "SoT structural loss (%)", ert_mean = "ERT mean (Ohm m)", ert_cv = "ERT CV", ert_pc1 = "ERT index (PC1)")
trees$sot_loss_pct <- trees$sot_structural_loss
kr_fit <- function(o, adjust) {
  sdx <- sd(unique(o[c("Tree", "x")])$x); if (!is.finite(sdx) || sdx == 0) return(c(b = NA, p = NA))
  o$x <- (o$x - mean(unique(o[c("Tree", "x")])$x)) / sdx
  f <- if (adjust) y ~ x + SPECIES + (1 | Tree) + (1 | date) else y ~ x + (1 | Tree) + (1 | date)
  m <- suppressMessages(lmer(f, data = o, REML = TRUE))
  c(b = fixef(m)[["x"]], p = anova(m, ddf = "Kenward-Roger")["x", "Pr(>F)"])
}
cell <- function(o, adjust = FALSE) {
  o <- o %>% filter(is.finite(x)); fit <- kr_fit(o, adjust)
  if (is.na(fit[["b"]])) return("–")
  loo <- sapply(unique(o$Tree), function(t) kr_fit(o %>% filter(Tree != t), adjust)[["b"]])
  sprintf("β = %.2f (%s)%s\nLOO %.2f to %.2f", fit[["b"]], ifelse(fit[["p"]] < 0.001, "<0.001", sprintf("%.3f", fit[["p"]])),
          ifelse(fit[["p"]] < 0.10, "*", ""), min(loo, na.rm = TRUE), max(loo, na.rm = TRUE))
}
obs <- function(d, m) fx %>% mutate(tree = as.character(Tree)) %>% inner_join(d %>% select(tree, x = all_of(m)), by = "tree")
grp <- list()
for (s in c("Wetland", "Upland")) {
  for (spp in c("bg", "rm", "hem", "ro")) {
    d <- trees %>% filter(site == s, SPECIES == spp); if (nrow(d) < 6) next
    grp[[length(grp) + 1]] <- c(group = sprintf("%s — %s", s, c(bg = "N. sylvatica", rm = "A. rubrum", hem = "T. canadensis", ro = "Q. rubra")[[spp]]),
      n = sprintf("%d / %d", nrow(d), sum(fx$Tree %in% d$tree)), sapply(names(metrics), function(m) cell(obs(d, m))))
  }
  d <- trees %>% filter(site == s)
  grp[[length(grp) + 1]] <- c(group = sprintf("%s — all species, species-adjusted", s), n = sprintf("%d / %d", nrow(d), sum(fx$site == s)),
    sapply(names(metrics), function(m) cell(obs(d, m), adjust = TRUE)))
}
si <- as.data.frame(do.call(rbind, grp), stringsAsFactors = FALSE)
names(si) <- c("Group", "Trees / measurements", unname(metrics))
write.csv(si, "outputs/tables/SI_decay_metric_correlations.csv", row.names = FALSE)
print(si)
