# ============================================================
# 07_manuscript_numbers.R
#
# Every number quoted in the manuscript text, computed in one place from
# the analysis dataset and the saved models, with the sample sizes,
# confidence intervals and standard errors requested by co-authors
# (docs/comment_tracker.md: Bradford #16/#17/#19/#27/#30/#31/#36/#40).
#
# Run after 4_analysis/ and 5_models/01-06.
#
# Outputs (outputs/tables/manuscript/):
#   manuscript_numbers.txt           human-readable summary for the text
#   site_summary.csv                 site means with n, SE, tree-bootstrap 95 % CI
#   species_tree_means.csv           tree-level species x site means, n trees, SE, 95 % CI
#   monthly_means.csv, summer_by_year.csv
#   variance_partitioning.csv        shares that sum to 100 %
#   coefficients_<model>.csv         SI tables: every fixed effect with SE, t, p
#   species_slopes_<model>.csv       species-specific standardized slopes with SE
#   vif_<model>.csv
#   model_data_info.csv              n obs, n trees, date range per model
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(lme4); library(car)
})
set.seed(20260930)
OUT <- file.path("outputs", "tables", "manuscript")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
txt <- character(0)
say <- function(...) { l <- paste0(...); txt <<- c(txt, l); cat(l, "\n") }
se <- function(x) sd(x) / sqrt(length(x))

# ------------------------------------------------------------
# Data
# ------------------------------------------------------------
d <- read.csv(file.path("data", "final", "stem_ch4_flux.csv"), stringsAsFactors = FALSE) %>%
  filter(!is.na(CH4_flux_nmolpm2ps), !is.na(PLOT)) %>%
  group_by(Tree) %>% mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>% ungroup() %>%
  mutate(location = ifelse(PLOT == "BGS", "Wetland", "Upland"), date = as.Date(date),
         species_full = dplyr::recode(SPECIES, bg = "N. sylvatica", hem = "T. canadensis", rm = "A. rubrum", ro = "Q. rubra"),
         month = as.integer(format(date, "%m")), year = as.integer(format(date, "%Y")),
         flux = CH4_flux_nmolpm2ps, asinh_flux = asinh(flux))   # flux already in nmol m-2 s-1: the model scale
dd <- sort(unique(d$date)); d$round <- cumsum(c(TRUE, diff(dd) > 7))[match(d$date, dd)]

say("=== DATA ===")
say(sprintf("Observations: %d; trees: %d; dates %s to %s; sampling rounds (>7 d apart): %d",
            nrow(d), n_distinct(d$Tree), min(d$date), max(d$date), max(d$round)))
say(sprintf("Flux source: %s", paste(names(table(d$flux_source)), table(d$flux_source), sep = " = ", collapse = ", ")))
say(sprintf("Analyzer: %s", paste(names(table(d$inst_label)), table(d$inst_label), sep = " = ", collapse = ", ")))

# tree-cluster bootstrap CI for a site mean (observations are nested in trees)
boot_ci <- function(df, B = 2000) {
  trees <- unique(df$Tree); byt <- split(df$flux, df$Tree)
  m <- replicate(B, mean(unlist(byt[sample(as.character(trees), replace = TRUE)])))
  quantile(m, c(0.025, 0.975))
}

# ------------------------------------------------------------
# Site summaries (#16, #17, #19)
# ------------------------------------------------------------
site <- d %>% group_by(location) %>% group_modify(~ {
  ci <- boot_ci(.x)
  tibble(n_obs = nrow(.x), n_trees = n_distinct(.x$Tree), n_rounds = n_distinct(.x$round),
         mean = mean(.x$flux), se = se(.x$flux), ci_lo = ci[1], ci_hi = ci[2],
         median = median(.x$flux), q25 = quantile(.x$flux, .25), q75 = quantile(.x$flux, .75),
         pct_positive = 100 * mean(.x$flux > 0), pct_negative = 100 * mean(.x$flux < 0),
         pct_below_MDF = 100 * mean(.x$CH4_below_MDF, na.rm = TRUE),
         pct_detected_negative = 100 * mean(.x$flux < 0 & !.x$CH4_below_MDF, na.rm = TRUE),
         variance = var(.x$flux))
}) %>% ungroup()
write_csv(site, file.path(OUT, "site_summary.csv"))
say("\n=== SITE MEANS (nmol m-2 s-1; SE over observations; 95 % CI = tree-cluster bootstrap) ===")
for (i in seq_len(nrow(site))) with(site[i, ], say(sprintf(
  "%s: mean %.2f ± %.2f (95%% CI %.2f–%.2f), median %.3f, IQR %.3f–%.3f; n = %d obs from %d trees in %d rounds; %.1f%% positive, %.1f%% negative; %.1f%% below MDF; %.1f%% detected negative",
  location, mean, se, ci_lo, ci_hi, median, q25, q75, n_obs, n_trees, n_rounds, pct_positive, pct_negative,
  pct_below_MDF, pct_detected_negative)))
w <- site %>% select(location, mean, variance) %>% pivot_wider(names_from = location, values_from = c(mean, variance))
say(sprintf("Wetland/upland mean ratio %.0f-fold; variance ratio %.0f-fold", w$mean_Wetland / w$mean_Upland,
            w$variance_Wetland / w$variance_Upland))

# ------------------------------------------------------------
# Tree-level species means (#16, #17)
# ------------------------------------------------------------
tm <- d %>% group_by(location, species_full, Tree) %>% summarise(n_obs = n(), mean_CH4 = mean(flux), .groups = "drop")
sp <- tm %>% group_by(location, species_full) %>%
  summarise(n_trees = n(), n_obs = sum(n_obs), mean = mean(mean_CH4), se = se(mean_CH4),
            ci_lo = mean - qt(.975, n_trees - 1) * se, ci_hi = mean + qt(.975, n_trees - 1) * se,
            median = median(mean_CH4), .groups = "drop")
write_csv(sp, file.path(OUT, "species_tree_means.csv"))
say("\n=== TREE-LEVEL MEANS BY SPECIES x SITE (mean ± SE of tree means; 95 % t-CI; n trees / n obs) ===")
for (i in seq_len(nrow(sp))) with(sp[i, ], say(sprintf("%s %s: %.3f ± %.3f (95%% CI %.3f–%.3f), median %.3f; %d trees, %d obs",
                                                       location, species_full, mean, se, ci_lo, ci_hi, median, n_trees, n_obs)))

# ------------------------------------------------------------
# Seasonality
# ------------------------------------------------------------
mm <- d %>% group_by(location, month) %>% summarise(n_obs = n(), n_trees = n_distinct(Tree), mean = mean(flux), se = se(flux), .groups = "drop")
write_csv(mm, file.path(OUT, "monthly_means.csv"))
say("\n=== MONTHLY PEAK / MINIMUM ===")
for (loc in c("Wetland", "Upland")) {
  x <- mm %>% filter(location == loc)
  hi <- x[which.max(x$mean), ]; lo <- x[which.min(x$mean), ]
  say(sprintf("%s: peak month %d = %.2f ± %.2f (n = %d); minimum month %d = %.3f (n = %d)",
              loc, hi$month, hi$mean, hi$se, hi$n_obs, lo$month, lo$mean, lo$n_obs))
}
sy <- d %>% filter(month %in% 6:8) %>% group_by(location, year) %>%
  summarise(n_obs = n(), mean = mean(flux), se = se(flux), .groups = "drop")
syb <- d %>% filter(month %in% 6:8, species_full == "N. sylvatica") %>% group_by(year) %>%
  summarise(n_obs = n(), mean = mean(flux), se = se(flux), .groups = "drop") %>% mutate(location = "N. sylvatica")
write_csv(bind_rows(sy, syb), file.path(OUT, "summer_by_year.csv"))
say("\n=== SUMMER (JUN–AUG) MEANS BY YEAR ===")
for (i in seq_len(nrow(sy))) with(sy[i, ], say(sprintf("%s %d: %.2f ± %.2f (n = %d)", location, year, mean, se, n_obs)))
for (i in seq_len(nrow(syb))) with(syb[i, ], say(sprintf("N. sylvatica %d: %.2f ± %.2f (n = %d)", year, mean, se, n_obs)))
rmx <- d %>% group_by(location, round) %>% summarise(date = min(date), mean = mean(flux), .groups = "drop") %>%
  group_by(location) %>% slice_max(mean, n = 1)
for (i in seq_len(nrow(rmx))) with(rmx[i, ], say(sprintf("Peak round mean %s: %.1f (%s)", location, mean, date)))

# ------------------------------------------------------------
# Upland detection by analyzer, with and without matching the 2025 conditions
# (LGR 2023-24 vs LI-7810 2025; matched = May-Oct and soil temperature (TS_Ha1) and shallow
# soil water (NEON) within the 10-90 % range of the 2025 measurements)
# ------------------------------------------------------------
env_u <- read.csv(file.path("data", "final", "environment_hourly.csv")) %>% transmute(hr = substr(datetime, 1, 13), TS_Ha1, NEON_SWC_shallow)
up <- read.csv(file.path("data", "final", "stem_ch4_flux.csv")) %>% filter(location == "Upland", !is.na(CH4_flux_nmolpm2ps)) %>%
  mutate(hr = sub(" ", "T", substr(sample_hour_est, 1, 13)), month = as.integer(substr(date, 6, 7))) %>% left_join(env_u, by = "hr")
rng <- up %>% filter(inst_label == "LI-7810") %>% summarise(t_lo = quantile(TS_Ha1, .1, na.rm = TRUE), t_hi = quantile(TS_Ha1, .9, na.rm = TRUE),
                                                          w_lo = quantile(NEON_SWC_shallow, .1, na.rm = TRUE), w_hi = quantile(NEON_SWC_shallow, .9, na.rm = TRUE))
up$matched <- up$month %in% 5:10 & between(up$TS_Ha1, rng$t_lo, rng$t_hi) & between(up$NEON_SWC_shallow, rng$w_lo, rng$w_hi)
det <- bind_rows(up %>% mutate(subset = "all"), up %>% filter(month %in% 5:10) %>% mutate(subset = "May-Oct"),
                 up %>% filter(matched %in% TRUE) %>% mutate(subset = "May-Oct, 2025 soil T and moisture range")) %>%
  group_by(subset, inst_label) %>% summarise(n = n(), pct_emission = 100 * mean(CH4_det_class == "emission"),
                                             pct_uptake = 100 * mean(CH4_det_class == "uptake"), median_MDF = median(CH4_MDF), .groups = "drop")
write_csv(det, file.path(OUT, "upland_detection_by_analyzer.csv"))
say("\n=== UPLAND DETECTION BY ANALYZER ===")
for (i in seq_len(nrow(det))) with(det[i, ], say(sprintf("%s, %s: n = %d, emission %.1f%%, uptake %.1f%%, median MDF %.3f", subset, inst_label, n, pct_emission, pct_uptake, median_MDF)))

# ------------------------------------------------------------
# Variance partitioning (#36, #37): shares of one model sum to 100 %
# ------------------------------------------------------------
say("\n=== VARIANCE PARTITIONING (asinh(flux, nmol m-2 s-1); random-intercept variance shares) ===")
vp <- list()
m0 <- lmer(asinh_flux ~ 1 + (1 | location), data = d)
v0 <- as.data.frame(VarCorr(m0)); tot0 <- sum(v0$vcov)
vp[[1]] <- tibble(model = "all obs: (1|site)", component = c("site", "residual"),
                  pct = 100 * c(v0$vcov[v0$grp == "location"], v0$vcov[v0$grp == "Residual"]) / tot0)
for (loc in c("Wetland", "Upland")) {
  m1 <- lmer(asinh_flux ~ 1 + (1 | species_full / Tree), data = d[d$location == loc, ])
  v <- as.data.frame(VarCorr(m1)); tot <- sum(v$vcov)
  vp[[loc]] <- tibble(model = paste(loc, ": (1|species/tree)"), component = c("species", "tree within species", "residual"),
                      pct = 100 * c(v$vcov[v$grp == "species_full"], v$vcov[grepl("Tree", v$grp)], v$vcov[v$grp == "Residual"]) / tot)
}
vp <- bind_rows(vp); write_csv(vp, file.path(OUT, "variance_partitioning.csv"))
for (mdl in unique(vp$model)) { x <- vp[vp$model == mdl, ]
  say(sprintf("%s — %s (sum %.0f%%)", mdl, paste(sprintf("%s %.1f%%", x$component, x$pct), collapse = ", "), sum(x$pct))) }

# ------------------------------------------------------------
# Models: SI coefficient tables, species slopes with SE, VIF, n and dates
# ------------------------------------------------------------
models <- list(
  wetland_core = "outputs/models/bgs_final/m_core_asinh.rds",
  wetland_full = "outputs/models/bgs_final/m_final.rds",
  upland_A = "outputs/models/ems_instantaneous/m_final.rds",
  upland_B = "outputs/models/ems_bgs_drivers/m_final.rds")
info_files <- c(wetland = "outputs/models/bgs_final/model_data_info.csv",
                upland_A = "outputs/models/ems_instantaneous/model_data_info.csv",
                upland_B = "outputs/models/ems_bgs_drivers/model_data_info.csv")
info <- bind_rows(lapply(names(info_files), function(k)
  if (file.exists(info_files[[k]])) read_csv(info_files[[k]], show_col_types = FALSE) %>% mutate(model = k, .before = 1)))
write_csv(info, file.path(OUT, "model_data_info.csv"))
say("\n=== MODEL DATA ===")
for (i in seq_len(nrow(info))) with(info[i, ], say(sprintf("%s: %d obs from %d trees, %s to %s", model, n_obs, n_trees,
                                                           as.Date(first), as.Date(last))))

# species-specific slope of each continuous term = base coefficient + species interaction,
# SE from the fixed-effect covariance matrix
species_slopes <- function(m) {
  mt <- as_lt(m)
  b <- fixef(m); V <- as.matrix(vcov(m)); nm <- names(b)
  comp <- strsplit(nm, ":", fixed = TRUE)
  is_sp <- function(x) grepl("^species", x)
  cont <- which(vapply(comp, function(x) !any(is_sp(x)) && !"(Intercept)" %in% x, logical(1)))
  sp_levels <- unique(unlist(lapply(comp, function(x) x[is_sp(x)])))
  mf <- model.frame(m); ref <- levels(factor(mf$species))[1]
  out <- list()
  for (k in cont) for (s in c(ref, sub("^species", "", sp_levels))) {
    w <- setNames(rep(0, length(b)), nm); w[k] <- 1
    if (s != ref) {
      j <- which(vapply(comp, function(x) setequal(x, c(comp[[k]], paste0("species", s))), logical(1)))
      if (length(j)) w[j] <- 1
    }
    est <- sum(w * b); s_e <- sqrt(as.numeric(t(w) %*% V %*% w))
    df <- if (!is.null(mt)) lmerTest::contest1D(mt, w)$df else Inf
    out[[length(out) + 1]] <- tibble(term = nm[k], species = s, estimate = est, se = s_e, df = df,
                                     t = est / s_e, p = 2 * pt(-abs(est / s_e), df))
  }
  bind_rows(out)
}
# p-values: Satterthwaite t-tests (lmerTest); falls back to normal approximation
as_lt <- function(m) tryCatch(lmerTest::as_lmerModLmerTest(m), error = function(e)
  tryCatch(suppressMessages(lmerTest::lmer(formula(m), data = m@frame, REML = isREML(m))), error = function(e2) NULL))
coef_table <- function(m) {
  mt <- as_lt(m)
  if (!is.null(mt)) { s <- summary(mt)$coefficients
    return(tibble(term = rownames(s), estimate = s[, 1], se = s[, 2], df = s[, "df"], t = s[, "t value"], p = s[, "Pr(>|t|)"])) }
  s <- summary(m)$coefficients
  tibble(term = rownames(s), estimate = s[, 1], se = s[, 2], df = Inf, t = s[, 3], p = 2 * pnorm(-abs(s[, 3])))
}
sp_names <- c(bg = "N. sylvatica", hem = "T. canadensis", rm = "A. rubrum", ro = "Q. rubra")
for (k in names(models)) {
  if (!file.exists(models[[k]])) { say("missing model: ", models[[k]]); next }
  m <- readRDS(models[[k]])
  ct <- coef_table(m); write_csv(ct, file.path(OUT, paste0("coefficients_", k, ".csv")))
  sl <- species_slopes(m) %>% mutate(species = dplyr::recode(species, !!!sp_names))
  write_csv(sl, file.path(OUT, paste0("species_slopes_", k, ".csv")))
  vv <- vif(m); vt <- tibble(term = rownames(vv), GVIF = vv[, "GVIF"], df = vv[, "Df"], GVIF_adj = vv[, "GVIF^(1/(2*Df))"])
  write_csv(vt, file.path(OUT, paste0("vif_", k, ".csv")))
  r2 <- as.numeric(performance::r2_nakagawa(m, tolerance = 1e-10)$R2_marginal)
  vc <- as.data.frame(VarCorr(m)); icc <- vc$vcov[vc$grp == "Tree"] / sum(vc$vcov)
  icc_date <- if (any(vc$grp == "date")) vc$vcov[vc$grp == "date"] / sum(vc$vcov) else NA
  say(sprintf("\n=== %s: R2 (Nakagawa marginal) = %.1f%%, AIC = %.1f, BIC = %.1f, k = %d, ICC = %.3f, n = %d ===",
              k, 100 * r2, AIC(m), BIC(m), length(fixef(m)), icc, nobs(m)))
  say(sprintf("  residual variance shares: tree %.1f%%, sampling date %.1f%%, residual %.1f%%", 100 * icc, 100 * icc_date,
              100 * vc$vcov[vc$grp == "Residual"] / sum(vc$vcov)))
  say(sprintf("  max GVIF = %.1f (%s); max GVIF^(1/(2Df)) = %.2f", max(vt$GVIF), vt$term[which.max(vt$GVIF)], max(vt$GVIF_adj)))
  for (i in seq_len(nrow(sl))) with(sl[i, ], say(sprintf("  %-45s %-14s %7.3f ± %.3f (p = %.2g)", term, species, estimate, se, p)))
}

# ------------------------------------------------------------
# Tomography (decay classes from the tomography paper; ERT-flux correlations)
# ------------------------------------------------------------
say("\n=== TOMOGRAPHY ===")
tc <- read_csv(file.path("data", "final", "tomography_classes.csv"), show_col_types = FALSE)
cls <- tc %>% filter(!is.na(decay_class))
say(sprintf("Decay classes (n = %d trees with SoT and ERT): %s", nrow(cls),
            paste(sprintf("%s %d (%.0f%%)", names(table(cls$decay_class)), table(cls$decay_class),
                          100 * prop.table(table(cls$decay_class))), collapse = "; ")))
bys <- round(100 * prop.table(table(cls$site, cls$decay_class), 1))
for (s in rownames(bys)) say(sprintf("  %s: %s", s, paste(colnames(bys), bys[s, ], "%", collapse = "; ")))
say(sprintf("SoT structural loss > 1%%: wetland %.0f%%, upland %.0f%%; mean ± SD wetland %.1f ± %.1f%%, upland %.1f ± %.1f%%",
            100 * mean(cls$sot_structural_loss[cls$site == "Wetland"] > 1), 100 * mean(cls$sot_structural_loss[cls$site == "Upland"] > 1),
            mean(cls$sot_structural_loss[cls$site == "Wetland"]), sd(cls$sot_structural_loss[cls$site == "Wetland"]),
            mean(cls$sot_structural_loss[cls$site == "Upland"]), sd(cls$sot_structural_loss[cls$site == "Upland"])))
ec <- file.path("outputs", "figures", "tomography", "ert_correlation_summary.csv")
if (file.exists(ec)) {
  e <- read_csv(ec, show_col_types = FALSE)
  for (met in c("ert_cv", "pc1_paper")) {
    x <- e %>% filter(metric == met)
    if (nrow(x)) say(sprintf("%s vs tree-mean CH4: %s", met,
                             paste(sprintf("%s r = %.2f (p = %.3g, n = %d)", x$group, x$r, x$p, x$n), collapse = "; ")))
  }
}
cvs <- read_csv(file.path("data", "package", "tomography", "tomography_results_compiled.csv"), show_col_types = FALSE) %>%
  mutate(tree = as.character(tree)) %>% inner_join(tc %>% mutate(tree = as.character(tree)), by = "tree")
for (s in c("Wetland", "Upland")) say(sprintf("ERT CV %s: %.3f ± %.3f (n = %d)", s,
  mean(cvs$ert_cv[cvs$site == s]), sd(cvs$ert_cv[cvs$site == s]), sum(cvs$site == s)))

# ------------------------------------------------------------
# 2025 drought (natural experiment): wetland water table and flux, August 2024 vs 2025
# ------------------------------------------------------------
say("\n=== 2025 DROUGHT ===")
al <- read.csv(file.path("data", "final", "environment_hourly.csv")) %>% mutate(dt = substr(datetime, 1, 13))
dr <- read.csv(file.path("data", "final", "stem_ch4_flux.csv")) %>% filter(location == "Wetland") %>%
  mutate(dt = substr(sub(" ", "T", sample_hour_est), 1, 13), yr = substr(date, 1, 4), m = as.integer(substr(date, 6, 7))) %>%
  left_join(al %>% select(dt, bvs_wtd_cm, bgs_wtd_cm), by = "dt")
dd <- dr %>% filter(m %in% c(8, 10), yr %in% c("2024", "2025")) %>% group_by(yr, m) %>%
  summarise(n = n(), wt_bvs = mean(bvs_wtd_cm, na.rm = TRUE), wt_bgs = mean(bgs_wtd_cm, na.rm = TRUE),
            flux = mean(CH4_flux_nmolpm2ps), flux_bg = mean(CH4_flux_nmolpm2ps[SPECIES == "bg"]), .groups = "drop")
write_csv(dd, file.path(OUT, "drought_2025.csv"))
for (i in seq_len(nrow(dd))) with(dd[i, ], say(sprintf("%s-%02d: n = %d; water table BVS %.0f cm, BGS %.0f cm; wetland mean %.2f, N. sylvatica %.2f nmol m-2 s-1",
                                                       yr, m, n, wt_bvs, wt_bgs, flux, flux_bg)))
up25 <- read.csv(file.path("data", "final", "stem_ch4_flux.csv")) %>% filter(location == "Upland", substr(date, 1, 4) == "2025") %>% pull(date) %>% unique() %>% sort()
say("Upland sampling dates 2025: ", paste(up25, collapse = ", "))

sy <- file.path("outputs", "tables", "synchrony_by_group.csv")   # from 4_analysis/04_repeatability.R
if (file.exists(sy)) {
  say("\n=== SYNCHRONY: date share of within-tree variance (asinh scale) ===")
  sy <- read.csv(sy)
  for (i in seq_len(nrow(sy))) with(sy[i, ], say(sprintf("%s %s: %.0f%% (LRT p = %.2g); tree %.0f%%, date %.0f%%, residual %.0f%% of total",
                                                       location, species_full, 100 * date_share_within_tree, date_lrt_p, 100 * tree_share, 100 * date_share, 100 * resid_share)))
}

mc <- file.path("outputs", "tables", "model_checks", "model_checks_summary.txt")
if (file.exists(mc)) { say("\n=== MODEL CHECKS (5_models/06_model_checks.R) ==="); for (l in readLines(mc)) say(l) }

writeLines(txt, file.path(OUT, "manuscript_numbers.txt"))
message("\nWrote ", file.path(OUT, "manuscript_numbers.txt"))
