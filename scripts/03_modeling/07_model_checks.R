# ============================================================
# 07_model_checks.R
#
# Diagnostics for the final driver models (01_bgs_model.R, 02_ems_model_A.R,
# 03_ems_model_B.R). Every check uses the saved final model and its data.
#
#   1. Variance explained: Nakagawa marginal/conditional R2 for species only,
#      core, and final model; how much the environmental drivers add beyond species
#   2. Sampling-date random effect: AIC with/without, SE inflation, p-values
#   3. Out-of-sample: leave-one-sampling-date-out and leave-one-tree-out
#      predictive R2; train 2023-24 / test 2025
#   4. Temperature vs season: share of the temperature predictor explained by the
#      seasonal cycle; does temperature explain flux beyond its seasonal cycle?
#   5. Response-transformation sensitivity of the species-specific terms
#   6. Upland Model A vs Model B on the observations both can use
#   7. Window screening null: date-block permutation of max |r| over all windows
#      (accounts for searching 112 windows and for shared sampling dates)
#
# Outputs: outputs/tables/model_checks/*.csv and model_checks_summary.txt
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(lme4); library(lmerTest); library(performance); library(RcppRoll)
})
set.seed(42)
OUT <- "outputs/tables/model_checks"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink_file <- file.path(OUT, "model_checks_summary.txt"); cat("", file = sink_file)
say <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = sink_file, append = TRUE) }

MODELS <- list(
  wetland = list(dir = "outputs/models/bgs_final",         site = "Wetland", label = "Wetland"),
  upland_A = list(dir = "outputs/models/ems_instantaneous", site = "Upland",  label = "Upland, Model A (upland drivers)"),
  upland_B = list(dir = "outputs/models/ems_bgs_drivers",   site = "Upland",  label = "Upland, Model B (wetland-type drivers)"))

r2n <- function(m) { r <- suppressWarnings(performance::r2_nakagawa(m, tolerance = 1e-10)); c(marg = as.numeric(r$R2_marginal), cond = as.numeric(r$R2_conditional)) }
pr2 <- function(o, p) { ok <- !is.na(p) & !is.na(o); 1 - sum((o[ok] - p[ok])^2) / sum((o[ok] - mean(o[ok]))^2) }
refit <- function(m, data, form = formula(m)) lmer(form, data = data, REML = FALSE, control = lmerControl(calc.derivs = FALSE))
preds_of <- function(m) setdiff(all.vars(formula(m)), c("CH4_flux_asinh", "species", "Tree", "date"))
core_form <- function(m) {
  p <- preds_of(m); ts <- p[grepl("^TS_|^s10t", p)][1]; mo <- setdiff(p[!grepl("^TS_|^s10t", p)], character())[1]
  as.formula(sprintf("CH4_flux_asinh ~ %s * %s * species + (1|Tree) + (1|date)", ts, mo))
}

fits <- list(); r2_rows <- list(); date_rows <- list(); oos_rows <- list(); seas_rows <- list(); tr_rows <- list()
for (k in names(MODELS)) {
  M <- MODELS[[k]]
  m <- readRDS(file.path(M$dir, "m_final.rds"))
  d <- readRDS(file.path(M$dir, "model_data_scaled.rds")) %>%
    mutate(Tree = factor(Tree), date = factor(as.Date(datetime)), species = factor(species))
  m <- refit(m, d); fits[[k]] <- list(m = m, d = d)
  p <- preds_of(m); ts <- p[grepl("^TS_|^s10t", p)][1]
  say("\n==================== ", M$label, " ====================")
  say("Final model: ", deparse1(formula(m)))
  say(sprintf("n = %d observations, %d trees, %d sampling dates, %s to %s", nobs(m), n_distinct(d$Tree),
              n_distinct(d$date), min(as.Date(d$datetime)), max(as.Date(d$datetime))))

  # ---- 1. variance explained ----
  m_sp <- refit(m, d, CH4_flux_asinh ~ species + (1|Tree) + (1|date))
  m_co <- refit(m, d, core_form(m))
  for (nm in c("species only", "core (temperature x moisture x species)", "final")) {
    mm <- switch(nm, "species only" = m_sp, "final" = m, m_co)
    r <- r2n(mm)
    r2_rows[[length(r2_rows) + 1]] <- tibble(model_set = k, model = nm, n = nobs(mm), k_fixed = length(fixef(mm)),
      R2_marginal = r[["marg"]], R2_conditional = r[["cond"]], AIC = AIC(mm), BIC = BIC(mm))
  }
  rr <- bind_rows(r2_rows) %>% filter(model_set == k)
  say(sprintf("Marginal R2: species only %.1f%%, core %.1f%%, final %.1f%% -> drivers add %.1f points beyond species",
              100 * rr$R2_marginal[1], 100 * rr$R2_marginal[2], 100 * rr$R2_marginal[3], 100 * (rr$R2_marginal[3] - rr$R2_marginal[1])))

  # ---- 2. sampling-date random effect ----
  f_nod <- update(formula(m), . ~ . - (1 | date))
  m_nod <- refit(m, d, f_nod)
  c1 <- summary(m)$coefficients; c0 <- summary(m_nod)$coefficients
  dr <- tibble(model_set = k, term = rownames(c1), estimate = c1[, "Estimate"], se = c1[, "Std. Error"], p = c1[, "Pr(>|t|)"],
               estimate_no_date = c0[rownames(c1), "Estimate"], se_no_date = c0[rownames(c1), "Std. Error"],
               p_no_date = c0[rownames(c1), "Pr(>|t|)"]) %>% mutate(se_ratio = se / se_no_date)
  date_rows[[k]] <- dr
  say(sprintf("Date random effect: AIC %.1f with vs %.1f without (dAIC = %.1f); SE ratio median %.2f; terms changing p<0.05 status: %s",
              AIC(m), AIC(m_nod), AIC(m) - AIC(m_nod), median(dr$se_ratio),
              paste(dr$term[(dr$p < 0.05) != (dr$p_no_date < 0.05)], collapse = ", ")))

  # ---- 3. out-of-sample ----
  y <- d$CH4_flux_asinh
  p_date <- rep(NA_real_, nrow(d)); p_date_sp <- p_date
  for (dt in levels(d$date)) {
    te <- d$date == dt; tr <- droplevels(d[!te, ])
    mt <- refit(m, tr); ms <- refit(m, tr, CH4_flux_asinh ~ species + (1|Tree) + (1|date))
    p_date[te] <- predict(mt, newdata = d[te, ], re.form = ~ (1 | Tree), allow.new.levels = TRUE)
    p_date_sp[te] <- predict(ms, newdata = d[te, ], re.form = ~ (1 | Tree), allow.new.levels = TRUE)
  }
  p_tree <- rep(NA_real_, nrow(d)); p_tree_sp <- p_tree
  for (tt in levels(d$Tree)) {
    te <- d$Tree == tt; tr <- droplevels(d[!te, ])
    mt <- refit(m, tr); ms <- refit(m, tr, CH4_flux_asinh ~ species + (1|Tree) + (1|date))
    p_tree[te] <- predict(mt, newdata = d[te, ], re.form = NA, allow.new.levels = TRUE)
    p_tree_sp[te] <- predict(ms, newdata = d[te, ], re.form = NA, allow.new.levels = TRUE)
  }
  yr <- as.numeric(format(as.Date(d$datetime), "%Y"))
  ts_split <- if (any(yr == 2025) && any(yr < 2025)) {
    tr <- droplevels(d[yr < 2025, ]); te <- d[yr == 2025, ]
    # restrict the test set to predictor values inside the training range (no extrapolation)
    inr <- Reduce(`&`, lapply(p, function(v) te[[v]] >= min(tr[[v]]) & te[[v]] <= max(tr[[v]])))
    mt <- refit(m, tr); ms <- refit(m, tr, CH4_flux_asinh ~ species + (1|Tree) + (1|date))
    c(n_test = nrow(te), n_in_range = sum(inr),
      all = pr2(te$CH4_flux_asinh, predict(mt, te, re.form = ~ (1 | Tree), allow.new.levels = TRUE)),
      in_range = pr2(te$CH4_flux_asinh[inr], predict(mt, te[inr, ], re.form = ~ (1 | Tree), allow.new.levels = TRUE)),
      species_in_range = pr2(te$CH4_flux_asinh[inr], predict(ms, te[inr, ], re.form = ~ (1 | Tree), allow.new.levels = TRUE)))
  } else c(n_test = 0, n_in_range = 0, all = NA, in_range = NA, species_in_range = NA)
  oos_rows[[k]] <- tibble(model_set = k,
    lodo_R2_final = pr2(y, p_date), lodo_R2_species = pr2(y, p_date_sp),
    loto_R2_final = pr2(y, p_tree), loto_R2_species = pr2(y, p_tree_sp),
    split2025_n = ts_split[["n_test"]], split2025_n_in_range = ts_split[["n_in_range"]],
    split2025_R2_all = ts_split[["all"]], split2025_R2_in_range = ts_split[["in_range"]],
    split2025_R2_species_in_range = ts_split[["species_in_range"]])
  o <- oos_rows[[k]]
  say(sprintf("Leave-one-date-out predictive R2: final %.2f vs species only %.2f; leave-one-tree-out: final %.2f vs species only %.2f",
              o$lodo_R2_final, o$lodo_R2_species, o$loto_R2_final, o$loto_R2_species))
  if (o$split2025_n > 0) say(sprintf("Train 2023-24 -> predict 2025: R2 %.2f on all %d obs; %.2f on the %d obs within the training range (species only %.2f)",
              o$split2025_R2_all, o$split2025_n, o$split2025_R2_in_range, o$split2025_n_in_range, o$split2025_R2_species_in_range))

  # ---- 4. temperature vs season ----
  if (!is.na(ts)) {
    doy <- as.numeric(format(as.Date(d$datetime), "%j"))
    d$s1 <- sin(2 * pi * doy / 365.25); d$c1 <- cos(2 * pi * doy / 365.25)
    sfit <- lm(as.formula(paste(ts, "~ s1 + c1")), data = d)
    d$ts_season <- fitted(sfit); d$ts_departure <- resid(sfit)
    f_seas <- as.formula(gsub(ts, "ts_season", deparse1(formula(m)), fixed = TRUE))
    m_seas <- refit(m, d, f_seas)
    m_dep  <- refit(m, d, update(f_seas, . ~ . + ts_departure * species))
    lrt <- anova(m_seas, m_dep)
    seas_rows[[k]] <- tibble(model_set = k, temperature_predictor = ts, R2_temperature_by_season = summary(sfit)$r.squared,
      AIC_final = AIC(m), AIC_seasonal_temperature = AIC(m_seas), AIC_plus_departure = AIC(m_dep),
      departure_LRT_chisq = lrt$Chisq[2], departure_LRT_df = lrt$Df[2], departure_LRT_p = lrt$`Pr(>Chisq)`[2])
    s <- seas_rows[[k]]
    say(sprintf("Temperature vs season: the seasonal cycle explains %.0f%% of %s. Model with seasonal temperature AIC %.1f vs final %.1f; adding the non-seasonal temperature departure (x species): chi2 = %.1f, df = %d, p = %.3g",
                100 * s$R2_temperature_by_season, ts, s$AIC_seasonal_temperature, s$AIC_final, s$departure_LRT_chisq, s$departure_LRT_df, s$departure_LRT_p))
  }

  # ---- 5. transformation sensitivity ----
  flux_nmol <- sinh(d$CH4_flux_asinh)
  tf <- list(`asinh(flux) [main]` = asinh(flux_nmol), `asinh(flux/0.1), more log-like` = asinh(flux_nmol / 0.1),
             `asinh(flux/10), more linear` = asinh(flux_nmol / 10),
             `rank-normal` = qnorm((rank(flux_nmol) - 0.5) / length(flux_nmol)))
  for (nm in names(tf)) {
    dd <- d; dd$CH4_flux_asinh <- tf[[nm]]
    cc <- summary(refit(m, dd))$coefficients
    tr_rows[[length(tr_rows) + 1]] <- tibble(model_set = k, transform = nm, term = rownames(cc),
                                             t = cc[, "t value"], p = cc[, "Pr(>|t|)"])
  }
}

r2_tab <- bind_rows(r2_rows); write_csv(r2_tab, file.path(OUT, "variance_explained.csv"))
write_csv(bind_rows(date_rows), file.path(OUT, "date_random_effect.csv"))
write_csv(bind_rows(oos_rows), file.path(OUT, "out_of_sample.csv"))
write_csv(bind_rows(seas_rows), file.path(OUT, "temperature_vs_season.csv"))
tr_tab <- bind_rows(tr_rows); write_csv(tr_tab, file.path(OUT, "transformation_sensitivity.csv"))
say("\nTransformation sensitivity (terms whose p < 0.05 status differs from the main asinh scale):")
tw <- tr_tab %>% mutate(sig = p < 0.05) %>% group_by(model_set, term) %>%
  summarise(main = sig[transform == "asinh(flux) [main]"], differs = paste(transform[sig != main], collapse = "; "), .groups = "drop") %>%
  filter(differs != "")
if (nrow(tw)) for (i in seq_len(nrow(tw))) say(sprintf("  %s  %s (main: %s) differs under: %s", tw$model_set[i], tw$term[i],
                                                       ifelse(tw$main[i], "p<0.05", "n.s."), tw$differs[i])) else say("  none")

# ---- 6. upland Model A vs Model B on shared observations ----
dA <- fits$upland_A$d; dB <- fits$upland_B$d
key <- intersect(paste(dA$Tree, dA$datetime), paste(dB$Tree, dB$datetime))
sA <- dA[paste(dA$Tree, dA$datetime) %in% key, ]; sB <- dB[paste(dB$Tree, dB$datetime) %in% key, ]
mA <- refit(fits$upland_A$m, droplevels(sA)); mB <- refit(fits$upland_B$m, droplevels(sB))
ab <- tibble(model = c("A (upland drivers)", "B (wetland-type drivers)"), n = c(nobs(mA), nobs(mB)),
             AIC = c(AIC(mA), AIC(mB)), BIC = c(BIC(mA), BIC(mB)), R2_marginal = c(r2n(mA)[["marg"]], r2n(mB)[["marg"]]))
write_csv(ab, file.path(OUT, "upland_A_vs_B_common_rows.csv"))
say(sprintf("\nUpland models on the %d observations both can use: A AIC %.1f, R2m %.1f%%; B AIC %.1f, R2m %.1f%%",
            length(key), ab$AIC[1], 100 * ab$R2_marginal[1], ab$AIC[2], 100 * ab$R2_marginal[2]))

# ---- 7. window screening null (date-block permutation of max |r| over windows) ----
a <- read_csv("data/processed/aligned_hourly_dataset.csv", show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC")) %>% arrange(datetime)
stopifnot(all(diff(as.numeric(a$datetime)) == 3600))
fl <- read_csv("data/processed/flux_with_quality_flags.csv", show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(format(as.POSIXct(sample_hour_est, tz = "UTC"), "%Y-%m-%d %H:%M:%S"), tz = "UTC"),
         date = as.Date(datetime), y = asinh(CH4_flux_nmolpm2ps)) %>%
  group_by(Tree) %>% mutate(y = y - mean(y)) %>% ungroup()        # as in 04_rolling_corrs.R
W <- seq(3, 14 * 24, by = 3)
rmw <- function(x, w) { m <- roll_mean(x, n = w, align = "right", fill = NA, na.rm = TRUE)
  ok <- roll_sum(as.numeric(!is.na(x)), n = w, align = "right", fill = NA) >= 0.5 * w; m[!ok %in% TRUE] <- NA; m }
perm_rows <- list()
for (k in names(MODELS)) {
  M <- MODELS[[k]]; pr <- preds_of(fits[[k]]$m); vars <- unique(sub("_raw_\\d+h$", "", pr[grepl("_raw_\\d+h$", pr)]))  # raw-mode drivers
  f <- fl %>% filter(location == M$site); idx <- match(f$datetime, a$datetime); ud <- unique(f$date)
  for (v in vars) {
    if (!v %in% names(a)) next
    X <- sapply(W, function(w) rmw(a[[v]], w))
    rr <- suppressWarnings(cor(f$y, X[idx, ], use = "pairwise.complete.obs"))[1, ]
    obs <- max(abs(rr), na.rm = TRUE); best <- W[which.max(abs(rr))]
    nulls <- replicate(500, {
      mp <- setNames(sample(ud), as.character(ud))
      sh <- as.numeric(difftime(mp[as.character(f$date)], f$date, units = "hours"))
      j <- idx + sh; j[j < 1 | j > nrow(X)] <- NA
      max(abs(suppressWarnings(cor(f$y, X[j, ], use = "pairwise.complete.obs"))), na.rm = TRUE)
    })
    near <- range(W[abs(rr) >= obs - 0.02], na.rm = TRUE)
    perm_rows[[length(perm_rows) + 1]] <- tibble(model_set = k, variable = v, best_window_h = best, max_abs_r = obs,
      r_range_over_windows = sprintf("%.2f to %.2f", min(rr, na.rm = TRUE), max(rr, na.rm = TRUE)),
      windows_within_0.02_h = sprintf("%d-%d", near[1], near[2]),
      null_median = median(nulls), null_95 = quantile(nulls, 0.95), p_perm = mean(nulls >= obs))
    say(sprintf("Window screening %s %s: best %d h, max|r| = %.2f; windows within 0.02 of the best: %d-%d h; date-permutation p = %.3f",
                k, v, best, obs, near[1], near[2], mean(nulls >= obs)))
  }
}
write_csv(bind_rows(perm_rows) %>% distinct(variable, best_window_h, .keep_all = TRUE), file.path(OUT, "window_permutation.csv"))
# ---- 8. wetland: core vs alternative and extended models on shared observations ----
# (the alternatives use drivers with different data coverage, e.g. tower latent heat ends
#  in Dec 2024, so they are compared on the observations every model can use)
dw <- fits$wetland$d; mw <- fits$wetland$m; pw <- preds_of(mw)
tsw <- pw[grepl("^TS_", pw)][1]; wtw <- pw[grepl("wtd", pw)][1]
pick <- function(pat) names(dw)[grepl(pat, names(dw))][1]
swc <- pick("^NEON_SWC_shallow_raw_"); le1 <- pick("^LE_Ha1_raw_")
ext_terms <- c(swc, pick("^FC_Ha1_anom_"), pick("^gcc_raw_"), pick("^tair_C_raw_"))
ext_terms <- ext_terms[!is.na(ext_terms)]
alt <- list(
  `Core: temperature x water table x species` = formula(mw),
  `Alternative: temperature x soil water content x species` = as.formula(sprintf("CH4_flux_asinh ~ %s * %s * species + (1|Tree) + (1|date)", tsw, swc)),
  `Alternative: latent heat x water table x species` = as.formula(sprintf("CH4_flux_asinh ~ %s * %s * species + (1|Tree) + (1|date)", le1, wtw)),
  `Season (day-of-year harmonics) x water table x species` = as.formula(sprintf("CH4_flux_asinh ~ (s1 + c1) * %s * species + (1|Tree) + (1|date)", wtw)),
  `Extended: core + selected additions x species` = update(formula(mw), as.formula(paste(". ~ . +", paste(paste0(ext_terms, " * species"), collapse = " + ")))))
doyw <- as.numeric(format(as.Date(dw$datetime), "%j")); dw$s1 <- sin(2 * pi * doyw / 365.25); dw$c1 <- cos(2 * pi * doyw / 365.25)
need <- unique(c(tsw, wtw, swc, le1, ext_terms))
dc <- droplevels(dw[complete.cases(dw[, need]), ])
cmp <- bind_rows(lapply(names(alt), function(nm) { mm <- refit(mw, dc, alt[[nm]])
  tibble(model = nm, n = nobs(mm), k_fixed = length(fixef(mm)), R2_marginal = r2n(mm)[["marg"]], AIC = AIC(mm), BIC = BIC(mm)) })) %>%
  mutate(dAIC = AIC - AIC[1], dBIC = BIC - BIC[1])
write_csv(cmp, file.path(OUT, "wetland_model_comparison_common_rows.csv"))
say(sprintf("\nWetland model comparison on the %d observations all models can use (%s to %s):", nrow(dc),
            min(as.Date(dc$datetime)), max(as.Date(dc$datetime))))
for (i in seq_len(nrow(cmp))) say(sprintf("  %-58s R2m %.1f%%  AIC %.1f (d %+.1f)  BIC %.1f (d %+.1f)  k = %d",
                                          cmp$model[i], 100 * cmp$R2_marginal[i], cmp$AIC[i], cmp$dAIC[i], cmp$BIC[i], cmp$dBIC[i], cmp$k_fixed[i]))

# species-specific slopes (reference coefficient + species interaction; Satterthwaite df)
sp_slopes <- function(mm, terms) {
  b <- fixef(mm); nm <- names(b); ref <- levels(mm@frame$species)[1]
  bind_rows(lapply(terms, function(tt) bind_rows(lapply(levels(mm@frame$species), function(sp) {
    w <- setNames(rep(0, length(b)), nm); w[tt] <- 1
    if (sp != ref) { ii <- paste0(tt, ":species", sp); alt_ii <- paste0("species", sp, ":", tt)
      if (ii %in% nm) w[ii] <- 1 else if (alt_ii %in% nm) w[alt_ii] <- 1
      # three-way terms are written temp:wtd:speciesX
    }
    ct <- lmerTest::contest1D(mm, w)
    tibble(term = tt, species = sp, estimate = ct$Estimate, se = ct$`Std. Error`, df = ct$df, p = ct$`Pr(>|t|)`) }))))
}
m_ext <- refit(mw, dc, alt[["Extended: core + selected additions x species"]])
m_cc  <- refit(mw, dc, alt[["Core: temperature x water table x species"]])
sl <- bind_rows(sp_slopes(m_cc, c(tsw, wtw, paste0(tsw, ":", wtw))) %>% mutate(model = "core (shared rows)"),
                sp_slopes(m_ext, c(tsw, wtw, paste0(tsw, ":", wtw))) %>% mutate(model = "extended (shared rows)"))
write_csv(sl, file.path(OUT, "wetland_species_slopes_core_vs_extended.csv"))
say("Species slopes, core vs extended (shared rows):")
for (i in seq_len(nrow(sl))) say(sprintf("  %-24s %-40s %-4s %6.3f ± %.3f (p = %.2g)", sl$model[i], sl$term[i], sl$species[i], sl$estimate[i], sl$se[i], sl$p[i]))

say("\nWritten to ", OUT)
