# ============================================================
# 01_tomography_classes.R
#
# Decay classes for the 60 stem-flux trees, reproduced exactly from the
# companion tomography paper (Thompson et al., "Internal decay in living
# trees: a quantitative tomography framework...";
# Tomography/Tree-Tomography/code/final_phase_and_scans.R):
#   - ERT axis: PC1 of 8 ERT metrics (mean, median, sd, cv, gini, entropy,
#     cma, radialgradient), each z-scored within species, PCA fitted on the
#     57 main-study trees with SoT; PC1 oriented so high = wetter/anomalous
#   - threshold: study-set mean of PC1 (anomalous moisture = PC1 > mean)
#   - SoT axis: structural loss (percent non-brown) > 1 %
#   - classes: I No Decay, II Incipient, III Active, IV Cavity
# SoT values are taken from the tomography study's tree table (Tree_ID_info.csv), which is
# upstream of tomography_results_compiled.csv (tree 433 differs: 16 % there vs 0 % in the
# compiled file). Both tables are in data/package/tomography/.
#
# Output: data/final/tomography_classes.csv
# ============================================================
suppressPackageStartupMessages(library(dplyr))
TOMO <- file.path("data", "package", "tomography")
tree_info <- read.csv(file.path(TOMO, "Tree_ID_info.csv"), stringsAsFactors = FALSE)
ert <- read.csv(file.path(TOMO, "ERT_application_results.csv"), stringsAsFactors = FALSE)
tree_info$tree <- as.character(tree_info$tree); ert$tree <- as.character(ert$tree)

m <- c("mean", "median", "sd", "cv", "gini", "entropy", "cma", "radialgradient")
train <- inner_join(tree_info, ert, by = "tree")          # 57 trees with SoT and ERT
stats <- train %>% group_by(species) %>%
  summarise(across(all_of(m), list(mu = ~ mean(.x, na.rm = TRUE), sg = ~ sd(.x, na.rm = TRUE))), .groups = "drop")
z <- function(d) {
  d <- left_join(d, stats, by = "species")
  X <- sapply(m, function(v) (d[[v]] - d[[paste0(v, "_mu")]]) / d[[paste0(v, "_sg")]])
  X[is.nan(X)] <- 0; X
}
pca <- prcomp(z(train), center = FALSE, scale. = FALSE)
flip <- if (pca$rotation["mean", 1] > 0) -1 else 1
thr <- mean(flip * (z(train) %*% pca$rotation)[, 1])

# all 60 flux trees with ERT (species from trees.csv for the 3 without SoT)
flux_sp <- read.csv(file.path("data", "package", "trees.csv"), stringsAsFactors = FALSE) %>%
  transmute(tree = as.character(Tree), species_flux = species, PLOT = plot)
allt <- ert %>% left_join(tree_info %>% select(tree, species, plot, percent_damaged, percent_solid_wood), by = "tree") %>%
  left_join(flux_sp, by = "tree") %>% mutate(species = coalesce(species, species_flux))
allt$ert_pc1 <- flip * (z(allt) %*% pca$rotation)[, 1]
out <- allt %>% transmute(
  tree, species, site = ifelse(coalesce(PLOT, plot) == "BGS", "Wetland", "Upland"),
  ert_pc1, ert_pc1_threshold = thr, sot_structural_loss = percent_damaged, sot_solid = percent_solid_wood,
  decay_class = case_when(is.na(percent_damaged) ~ NA_character_,
                          percent_damaged <= 1 & ert_pc1 <= thr ~ "I: No Decay",
                          percent_damaged <= 1 & ert_pc1 >  thr ~ "II: Incipient",
                          percent_damaged >  1 & ert_pc1 >  thr ~ "III: Active",
                          TRUE ~ "IV: Cavity"))
message("PC1 variance explained: ", round(100 * summary(pca)$importance[2, 1], 1), " %; threshold ", signif(thr, 3))
print(round(flip * pca$rotation[, 1], 2))
print(table(out$decay_class, useNA = "ifany"))
print(round(100 * prop.table(table(out$site, out$decay_class), 1)))
write.csv(out, file.path("data", "final", "tomography_classes.csv"), row.names = FALSE)
message("Wrote data/final/tomography_classes.csv (", nrow(out), " trees)")
