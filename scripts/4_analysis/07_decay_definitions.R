# ============================================================
# 07_decay_definitions.R
#
# Table S2: does the wood-condition / CH4 flux association depend on the metric?
# Three tree-level metrics against every CH4 measurement:
#   SoT structural loss  % of the SoT cross-section in non-brown (low-velocity) classes (PiCUS Q74)
#   ERT CV               coefficient of variation of resistivity (heterogeneity of moisture)
#   ERT index (PC1)      first principal component of eight ERT metrics, each z-scored within
#                        species (tomography paper; higher = wetter, more heterogeneous)
# Each metric is measured once per tree, so the tree is the unit of replication: every
# measurement enters flux ~ metric (per SD) + (1 | tree) + (1 | date), tested with
# Kenward-Roger degrees of freedom (about 8 for 10 trees), as in Figure 5 (06_tomography_flux.R).
# Rows: each species x site (10 trees), and per site all species with species as a fixed effect.
# Cells: change in flux (nmol m-2 s-1) per SD of the metric (p), and its range when each tree is
# omitted in turn. A last column repeats the ERT CV test on asinh(flux / 0.1), the scale of the
# environmental-driver models, which weights typical rather than mean flux.
#
# Inputs : data/package/tomography/ERT_application_results.csv, data/final/tomography_classes.csv,
#          data/final/stem_ch4_flux.csv
# Output : outputs/tables/SI_decay_metric_correlations.csv (read by 6_manuscript/01_build_si.R)
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(lmerTest) })

# ERT metrics per tomogram (tomography dataset), with the metric names used below
read_ert <- function(f) {
  x <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  m <- c("mean", "median", "sd", "cv", "gini", "entropy", "cma", "radialgradient")
  stopifnot(all(m %in% names(x)))
  names(x)[match(m, names(x))] <- paste0("ert_", sub("radialgradient", "radialgradiant", m))
  x
}
tc <- read.csv("data/final/tomography_classes.csv", stringsAsFactors = FALSE)
ert <- read_ert("data/package/tomography/ERT_application_results.csv")
fx <- read.csv("data/final/stem_ch4_flux.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(PLOT), !is.na(CH4_flux_nmolpm2ps)) %>%
  group_by(Tree) %>% mutate(SPECIES = first(na.omit(SPECIES))) %>% ungroup() %>%
  mutate(site = ifelse(PLOT == "BGS", "Wetland", "Upland"), y = CH4_flux_nmolpm2ps)

# one row per tree: ERT CV, ERT PC1 and SoT loss (tomography_classes.csv), species and site
trees <- ert %>% transmute(tree = as.character(tree), ert_cv) %>%
  inner_join(tc %>% mutate(tree = as.character(tree)), by = "tree") %>%
  inner_join(fx %>% group_by(tree = as.character(Tree)) %>%
               summarise(flux = mean(CH4_flux_nmolpm2ps), SPECIES = first(SPECIES), .groups = "drop"), by = "tree")

metrics <- c(sot_loss_pct = "SoT structural loss (%)", ert_cv = "ERT CV", ert_pc1 = "ERT index (PC1)")
trees$sot_loss_pct <- trees$sot_structural_loss
kr_fit <- function(o, adjust, test = TRUE) {   # test = FALSE: slope only (leave-one-out refits)
  sdx <- sd(unique(o[c("Tree", "x")])$x); if (!is.finite(sdx) || sdx == 0) return(c(b = NA, p = NA))
  o$x <- (o$x - mean(unique(o[c("Tree", "x")])$x)) / sdx
  f <- if (adjust) y ~ x + SPECIES + (1 | Tree) + (1 | date) else y ~ x + (1 | Tree) + (1 | date)
  m <- suppressMessages(lmer(f, data = o, REML = TRUE))
  c(b = fixef(m)[["x"]], p = if (test) anova(m, ddf = "Kenward-Roger")["x", "Pr(>F)"] else NA)
}
cell <- function(o, adjust = FALSE) {
  o <- o %>% filter(is.finite(x)); fit <- kr_fit(o, adjust)
  if (is.na(fit[["b"]])) return("–")
  loo <- sapply(unique(o$Tree), function(t) kr_fit(o %>% filter(Tree != t), adjust, test = FALSE)[["b"]])
  g <- function(v) formatC(v, digits = 2, format = "fg", flag = "#")
  sprintf("%s (%s)%s\nLOO %s to %s", g(fit[["b"]]), ifelse(fit[["p"]] < 0.001, "<0.001", sprintf("%.3f", fit[["p"]])),
          ifelse(fit[["p"]] < 0.10, "*", ""), g(min(loo, na.rm = TRUE)), g(max(loo, na.rm = TRUE)))
}
obs <- function(d, m) fx %>% mutate(tree = as.character(Tree)) %>% inner_join(d %>% select(tree, x = all_of(m)), by = "tree")
cell_asinh <- function(o, adjust = FALSE) {   # same test, response asinh(flux / 0.1)
  o <- o %>% filter(is.finite(x)) %>% mutate(y = asinh(y / 0.1)); fit <- kr_fit(o, adjust)
  if (is.na(fit[["b"]])) return("–")
  sprintf("%s (%s)%s", formatC(fit[["b"]], digits = 2, format = "fg", flag = "#"),
          ifelse(fit[["p"]] < 0.001, "<0.001", sprintf("%.3f", fit[["p"]])), ifelse(fit[["p"]] < 0.10, "*", ""))
}
grp <- list()
for (s in c("Wetland", "Upland")) {
  for (spp in c("bg", "rm", "hem", "ro")) {
    d <- trees %>% filter(site == s, SPECIES == spp); if (nrow(d) < 6) next
    grp[[length(grp) + 1]] <- c(group = sprintf("%s — %s", s, c(bg = "N. sylvatica", rm = "A. rubrum", hem = "T. canadensis", ro = "Q. rubra")[[spp]]),
      n = sprintf("%d / %d", nrow(d), sum(fx$Tree %in% d$tree)), sapply(names(metrics), function(m) cell(obs(d, m))),
      asinh = cell_asinh(obs(d, "ert_cv")))
  }
  d <- trees %>% filter(site == s)
  grp[[length(grp) + 1]] <- c(group = sprintf("%s — all species, species-adjusted", s), n = sprintf("%d / %d", nrow(d), sum(fx$site == s)),
    sapply(names(metrics), function(m) cell(obs(d, m), adjust = TRUE)), asinh = cell_asinh(obs(d, "ert_cv"), adjust = TRUE))
}
si <- as.data.frame(do.call(rbind, grp), stringsAsFactors = FALSE)
names(si) <- c("Group", "Trees / measurements", unname(metrics), "ERT CV, asinh(flux / 0.1)")
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(si, "outputs/tables/SI_decay_metric_correlations.csv", row.names = FALSE)
print(si)
