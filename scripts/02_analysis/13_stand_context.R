# ============================================================
# 13_stand_context.R
# Stand context for Table 1: share of live basal area and stand DBH quartiles of each focal
# species, total live basal area, and the DBH of the study trees.
#   Wetland (Black Gum Swamp): June 2025 variable-radius (prism) survey, 40 plots inside the
#     swamp outline; basal-area factor 10 ft2/acre = 2.2957 m2/ha per tallied stem (from the
#     survey workbook, Matthes_Lab/bgs-VRP-sampling/VRP-sampling.xlsx). Basal-area share =
#     tally share; DBH quartiles weight each stem by its expansion factor (proportional to 1/DBH^2).
#     Shares are also given for the plots within 75 m of the study trees.
#   Upland (EMS): 2019 Harvard Forest ForestGEO census (HF253 v6; Orwig, Foster & Ellison 2024),
#     live stems, whole 35-ha plot.
# Method follows the companion tomography study (Tree-Tomography/analysis/revision/scripts/stand_context.R).
# Outputs: outputs/tables/stand_context.csv, outputs/tables/table1_study_trees.csv,
#          outputs/tables/stand_totals.csv
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(sf) })
S <- "data/input/spatial"
BAF <- 10 * 0.09290304 / 0.40468564   # 10 ft2/acre in m2/ha
latin <- c(acerru = "A. rubrum", nysssy = "N. sylvatica", querru = "Q. rubra", tsugca = "T. canadensis")
code <- c(rm = "acerru", bg = "nysssy", ro = "querru", hem = "tsugca")

vrp_all <- read.csv(file.path(S, "BGS_VRP_2025.csv"))
vrp <- vrp_all %>% filter(status != "dead", !is.na(dbh))                 # alive + stressed
cen <- read.csv(file.path(S, "hf253-06-stems-2019.csv")) %>% filter(status == "A", !is.na(dbh))   # DBH in cm
n_plots <- length(unique(vrp_all$plot))

# plots within 75 m of the study trees
pl <- vrp_all %>% distinct(plot, lat, long) %>% st_as_sf(coords = c("long", "lat"), crs = 4326) %>% st_transform(26986)
tr <- st_read(file.path(S, "BGS_editable.geojson"), quiet = TRUE) %>% st_zm() %>% st_transform(26986)
near <- unique(pl$plot[apply(st_distance(pl, tr), 1, min) <= 75])

ems_ba <- cen %>% mutate(ba = pi * (dbh / 200)^2)            # m2 per stem
plot_area_ha <- 35
totals <- data.frame(
  site = c("Wetland", "Upland"),
  live_basal_area_m2_ha = c(nrow(vrp) * BAF / n_plots, sum(ems_ba$ba) / plot_area_ha),
  source = c(sprintf("Prism survey, %d plots, BAF 10 ft2/acre", n_plots), "ForestGEO 2019 census, 35 ha"))
share_bgs <- table(vrp$species) / nrow(vrp)
share_bgs_near <- { x <- vrp %>% filter(plot %in% near); table(x$species) / nrow(x) }
share_ems <- tapply(ems_ba$ba, ems_ba$sp, sum) / sum(ems_ba$ba)
wq <- function(x, w, p) { o <- order(x); x <- x[o]; cw <- cumsum(w[o]) / sum(w); sapply(p, function(pp) x[which(cw >= pp)[1]]) }

fx <- read.csv("data/processed/flux_with_quality_flags.csv")
ours <- fx %>% group_by(Tree) %>% summarise(site = first(location), sp = code[first(na.omit(SPECIES))], dbh = first(na.omit(DBH)), .groups = "drop")

rows <- list()
for (site in c("Wetland", "Upland")) for (s in names(latin)) {
  o <- ours$dbh[ours$sp == s & ours$site == site]; if (!length(o)) next
  if (site == "Wetland") { x <- vrp$dbh[vrp$species == s & vrp$dbh >= 10]; w <- 1 / x^2; sh <- share_bgs[[s]]; shn <- share_bgs_near[[s]] }
  else { x <- cen$dbh[cen$sp == s & cen$dbh >= 10]; w <- rep(1, length(x)); sh <- share_ems[[s]]; shn <- NA }
  q <- wq(x, w, c(.25, .5, .75))
  rows[[length(rows) + 1]] <- data.frame(site = site, species = latin[[s]], n_trees = length(o),
    dbh_mean_cm = round(mean(o), 1), dbh_min_cm = round(min(o), 1), dbh_max_cm = round(max(o), 1),
    pct_live_basal_area = round(100 * sh), pct_live_basal_area_near_trees = round(100 * shn),
    stand_dbh_p25_cm = round(q[1]), stand_dbh_p50_cm = round(q[2]), stand_dbh_p75_cm = round(q[3]))
}
tab <- bind_rows(rows)
write.csv(tab, "outputs/tables/stand_context.csv", row.names = FALSE)
write.csv(totals, "outputs/tables/stand_totals.csv", row.names = FALSE)
# dominant species (all species) for the site description
top <- function(sh, k = 5) { sh <- sort(sh, decreasing = TRUE)[1:k]; paste(sprintf("%s %d%%", names(sh), round(100 * sh)), collapse = ", ") }
writeLines(c(sprintf("Wetland live BA %.1f m2/ha (%d plots); all plots: %s; plots within 75 m of study trees (%d): %s",
                     totals$live_basal_area_m2_ha[1], n_plots, top(share_bgs), length(near), top(share_bgs_near)),
             sprintf("Upland live BA %.1f m2/ha (ForestGEO 35 ha): %s", totals$live_basal_area_m2_ha[2], top(share_ems))),
           "outputs/tables/stand_composition.txt")
print(tab); cat(readLines("outputs/tables/stand_composition.txt"), sep = "\n")
