# ============================================================
# 10_decay_binned.R  (exploratory)
#
# Binned and threshold versions of the wood-condition vs CH4 flux analysis, as an
# alternative to the continuous correlations of 06/07. Unit of analysis = tree
# (tree-mean asinh flux), as in the continuous analysis. For each species x site and
# for each site pooled (species-centred), and each metric:
#   median split   high vs low half: difference in mean, Wilcoxon p
#   tertiles       Kruskal-Wallis p across thirds; Spearman of tertile vs flux
#   class          tomography-study classes (decay II-IV vs I; moisture anomaly), Wilcoxon p
#   best split     threshold maximizing |t| between groups (>= 3 trees per side);
#                  p from 2,000 permutations of the metric among trees, repeating the
#                  whole threshold search each time (corrects for choosing the cut)
#   hinge          flux = a + b * max(0, metric - c), best c by least squares;
#                  permutation p as above
# Output: outputs/tables/decay_binned.csv
# ============================================================
suppressPackageStartupMessages({ library(dplyr) })
set.seed(42)
NPERM <- 2000

tc  <- read.csv("data/final/tomography_classes.csv", stringsAsFactors = FALSE)
ert <- read.csv("data/package/tomography/tomography_results_compiled.csv", stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
fx  <- read.csv("data/final/stem_ch4_flux.csv", stringsAsFactors = FALSE)
zs  <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
trees <- ert %>% transmute(tree = as.character(tree), ert_cv, ert_mean, ert_cma) %>%
  inner_join(tc %>% mutate(tree = as.character(tree)), by = "tree") %>%
  inner_join(fx %>% group_by(tree = as.character(Tree)) %>%
               summarise(y = mean(asinh(CH4_flux_nmolpm2ps)), SPECIES = first(SPECIES), site = first(location), .groups = "drop"),
             by = c("tree", "site")) %>%
  mutate(ert_cv_z = zs(ert_cv, SPECIES), ert_mean_z = zs(ert_mean, SPECIES), sot_loss = sqrt(sot_structural_loss),
         decay_any = decay_class != "I: No Decay", moisture_anom = ert_pc1 > ert_pc1_threshold)
METRICS <- c("ert_cv", "ert_cv_z", "ert_pc1", "sot_loss", "ert_mean_z", "ert_cma")

best_split <- function(x, y, min_n = 3) {
  o <- order(x); x <- x[o]; y <- y[o]; n <- length(x); best <- c(stat = 0, cut = NA, diff = NA)
  for (k in min_n:(n - min_n)) {
    if (x[k] == x[k + 1]) next
    a <- y[1:k]; b <- y[(k + 1):n]
    t <- (mean(b) - mean(a)) / sqrt(var(a) / length(a) + var(b) / length(b))
    if (is.finite(t) && abs(t) > abs(best[["stat"]])) best <- c(stat = t, cut = (x[k] + x[k + 1]) / 2, diff = mean(b) - mean(a))
  }
  best
}
hinge <- function(x, y, min_n = 3) {
  cs <- sort(x)[min_n:(length(x) - min_n)]; best <- c(rss = Inf, cut = NA, slope = NA)
  for (c0 in cs) { h <- pmax(0, x - c0); f <- lm(y ~ h); r <- sum(resid(f)^2)
    if (r < best[["rss"]]) best <- c(rss = r, cut = c0, slope = coef(f)[[2]]) }
  best[["r2"]] <- 1 - best[["rss"]] / sum((y - mean(y))^2); best
}
perm_p <- function(x, y, fun, stat) {
  obs <- stat(fun(x, y)); null <- replicate(NPERM, stat(fun(sample(x), y)))
  (1 + sum(null >= obs - 1e-12)) / (NPERM + 1)
}
one <- function(d, label) {
  out <- list()
  for (m in METRICS) {
    e <- d[!is.na(d[[m]]), ]; x <- e[[m]]; y <- e$y; n <- nrow(e); if (n < 8) next
    hi <- x > median(x)
    ter <- cut(rank(x, ties.method = "first"), 3, labels = FALSE)
    bs <- best_split(x, y); hg <- hinge(x, y)
    out[[m]] <- tibble(group = label, metric = m, n_trees = n,
      r = cor(x, y), r_p = cor.test(x, y)$p.value,
      median_diff = mean(y[hi]) - mean(y[!hi]), median_wilcox_p = suppressWarnings(wilcox.test(y[hi], y[!hi])$p.value),
      tertile_kw_p = kruskal.test(y, ter)$p.value, tertile_rho = cor(ter, y, method = "spearman"),
      split_cut = bs[["cut"]], split_diff = bs[["diff"]], split_n_high = sum(x > bs[["cut"]]),
      split_perm_p = perm_p(x, y, best_split, function(b) abs(b[["stat"]])),
      hinge_cut = hg[["cut"]], hinge_slope = hg[["slope"]], hinge_r2 = hg[["r2"]],
      hinge_perm_p = perm_p(x, y, hinge, function(b) b[["r2"]]))
  }
  for (cl in c("decay_any", "moisture_anom")) {
    g <- d[[cl]]; if (sum(g, na.rm = TRUE) < 2 || sum(!g, na.rm = TRUE) < 2) next
    out[[cl]] <- tibble(group = label, metric = cl, n_trees = sum(!is.na(g)), n_high = sum(g, na.rm = TRUE),
      median_diff = mean(d$y[g %in% TRUE]) - mean(d$y[g %in% FALSE]),
      median_wilcox_p = suppressWarnings(wilcox.test(d$y[g %in% TRUE], d$y[g %in% FALSE])$p.value))
  }
  bind_rows(out)
}
res <- list()
for (s in c("Wetland", "Upland")) {
  for (sp in unique(trees$SPECIES[trees$site == s])) res[[paste(s, sp)]] <- one(trees %>% filter(site == s, SPECIES == sp), paste(s, sp))
  pooled <- trees %>% filter(site == s) %>% group_by(SPECIES) %>% mutate(y = y - mean(y)) %>% ungroup()   # species-centred
  res[[paste(s, "pooled")]] <- one(pooled, paste(s, "pooled (species-centred)"))
}
res <- bind_rows(res)
write.csv(res, "outputs/tables/decay_binned.csv", row.names = FALSE)
options(width = 220)
print(as.data.frame(res %>% select(group, metric, n_trees, r, r_p, median_diff, median_wilcox_p, tertile_kw_p, split_cut, split_n_high, split_perm_p, hinge_r2, hinge_perm_p) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 2)))))
