# ============================================================
# 11_site_map.R
# Study-site map: (a) location in New England, (b) Prospect Hill terrain with the two
# sites and the flux towers, (c) wetland trees at Black Gum Swamp, (d) upland trees
# near the EMS tower. Tree positions are field GPS surveys (BGS 2025-04-17, EMS 2026-01-16).
#
# Inputs (data/input/spatial/, see data/README.md for sources):
#   BGS_editable.geojson, EMS_trees_editable.geojson   study-tree positions (WGS84)
#   elevation_ned/ (30-m NED DEM, MA State Plane), elevation_contours.geojson
#   tracts.shp (Harvard Forest tracts), Black_Gum_Swamp.kmz (swamp outline)
#   cb_2018_us_state_500k.zip (US Census state outlines)
# Tower coordinates: AmeriFlux site metadata (US-Ha1 EMS, US-Ha2 hemlock, US-xHA NEON).
# Output: outputs/figures/site_map.png / .pdf; outputs/tables/study_tree_coordinates.csv
# ============================================================
suppressPackageStartupMessages({
  library(sf); library(terra); library(ggplot2); library(dplyr); library(patchwork); library(ggspatial)
})
S <- "data/input/spatial"
OUT <- "outputs/figures"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CRS_MA <- 26986   # NAD83 / Massachusetts Mainland (m)

species_colors <- c("N. sylvatica" = "#2A7F7A", "Q. rubra" = "#6E8B3D", "A. rubrum" = "#A7DAD1", "T. canadensis" = "#C9D6A4")
spp_shapes <- c("N. sylvatica" = 22, "A. rubrum" = 24, "T. canadensis" = 21, "Q. rubra" = 23)
site_colors <- c("Wetland (BGS)" = "#2A7F7A", "Upland (EMS)" = "#6E8B3D")
sp_name <- c(BG = "N. sylvatica", RM = "A. rubrum", HEM = "T. canadensis", RO = "Q. rubra")

# ---- data ----
trees <- bind_rows(
  st_read(file.path(S, "BGS_editable.geojson"), quiet = TRUE) %>% mutate(site = "Wetland (BGS)"),
  st_read(file.path(S, "EMS_trees_editable.geojson"), quiet = TRUE) %>% mutate(site = "Upland (EMS)")) %>%
  st_zm() %>% transmute(tree = Name, species = sp_name[toupper(species)], site,
                        microsite = ifelse(site == "Wetland (BGS)", description, NA)) %>% st_transform(CRS_MA)
stopifnot(nrow(trees) == 60, !anyNA(trees$species))
write.csv(cbind(st_drop_geometry(trees), st_coordinates(st_transform(trees, 4326))) %>% rename(lon = X, lat = Y),
          "outputs/tables/study_tree_coordinates.csv", row.names = FALSE)

towers <- st_as_sf(data.frame(
  name = c("EMS tower (US-Ha1)", "Hemlock tower (US-Ha2)", "NEON tower (US-xHA)"),
  lon = c(-72.1715, -72.1779, -72.17265), lat = c(42.5378, 42.5393, 42.53691)),
  coords = c("lon", "lat"), crs = 4326) %>% st_transform(CRS_MA)

kdir <- file.path(tempdir(), "bgs_kmz"); unzip(file.path(S, "Black_Gum_Swamp.kmz"), exdir = kdir)
swamp <- st_read(list.files(kdir, pattern = "\\.kml$", full.names = TRUE)[1], quiet = TRUE) %>% st_zm() %>%
  st_transform(CRS_MA) %>% st_union()
swamp_poly <- tryCatch(st_cast(st_polygonize(swamp), "POLYGON"), error = function(e) NULL)

tracts <- st_read(file.path(S, "tracts.shp"), quiet = TRUE) %>% st_transform(CRS_MA)
dem <- rast(file.path(S, "elevation_ned"))
crs(dem) <- paste0("EPSG:", CRS_MA)

# ---- panel b extent: both sites and the towers, padded ----
bb <- st_bbox(st_union(st_geometry(trees), st_geometry(towers)))
pad <- 450
ext_b <- ext(bb["xmin"] - pad, bb["xmax"] + pad, bb["ymin"] - pad, bb["ymax"] + pad)
eb <- unname(as.vector(ext_b))   # xmin, xmax, ymin, ymax
dem_b <- crop(dem, ext_b)
hs <- shade(terrain(dem_b, "slope", unit = "radians"), terrain(dem_b, "aspect", unit = "radians"), 40, 315)
hs_df <- as.data.frame(hs, xy = TRUE); names(hs_df)[3] <- "hs"
el_df <- as.data.frame(dem_b, xy = TRUE); names(el_df)[3] <- "elev"
cont <- st_read(file.path(S, "elevation_contours.geojson"), quiet = TRUE) %>% st_transform(CRS_MA) %>%
  st_crop(st_bbox(c(xmin = eb[1], xmax = eb[2], ymin = eb[3], ymax = eb[4]), crs = st_crs(CRS_MA)))
site_hull <- trees %>% group_by(site) %>% summarise(.groups = "drop") %>% st_convex_hull() %>% st_buffer(25)

base_theme <- theme_bw(base_size = 10) +
  theme(axis.title = element_blank(), panel.grid = element_blank(), axis.text = element_text(size = 7),
        plot.tag = element_text(face = "bold", size = 12), plot.title = element_text(size = 10, face = "bold"))

pb <- ggplot() +
  geom_raster(data = el_df, aes(x, y, fill = elev), alpha = 0.9) +
  scale_fill_gradientn(colours = c("#F2F4EE", "#DCE3D0", "#BFCBA8", "#9FAE86"), name = "Elevation (m)") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = hs_df, aes(x, y, alpha = hs), fill = "grey15", show.legend = FALSE) +
  scale_alpha_continuous(range = c(0.35, 0)) +
  geom_sf(data = cont, colour = "grey45", linewidth = 0.15, alpha = 0.6) +
  { if (!is.null(swamp_poly)) geom_sf(data = swamp_poly, fill = "#2A7F7A", alpha = 0.25, colour = NA) } +
  geom_sf(data = swamp, colour = "#1F5E5A", linewidth = 0.5) +
  geom_sf(data = trees, aes(fill = site), shape = 21, colour = "white", stroke = 0.2, size = 1.8) +
  geom_sf(data = towers, shape = 24, fill = "black", colour = "white", size = 2.6) +
  geom_sf_text(data = towers, aes(label = sub(" \\(", "\n(", name)), size = 2.3, nudge_y = c(-80, 85, -95), nudge_x = c(150, 0, -40), lineheight = 0.9) +
  scale_fill_manual(values = site_colors, name = NULL, aesthetics = "fill") +
  coord_sf(xlim = eb[1:2], ylim = eb[3:4], crs = st_crs(CRS_MA), datum = st_crs(4326), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.6) +
  annotation_north_arrow(location = "tr", height = unit(0.8, "cm"), width = unit(0.6, "cm"), style = north_arrow_minimal()) +
  labs(tag = "b", title = "Prospect Hill tract, Harvard Forest") + base_theme +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "cm"), legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))

zoom <- function(site_name, tag, title) {
  t <- trees %>% filter(site == site_name)
  b <- st_bbox(st_buffer(st_union(t), 30))
  g <- ggplot()
  if (site_name == "Wetland (BGS)") {
    if (!is.null(swamp_poly)) g <- g + geom_sf(data = swamp_poly, fill = "#2A7F7A", alpha = 0.12, colour = NA)
    g <- g + geom_sf(data = swamp, colour = "#1F5E5A", linewidth = 0.4)
  }
  zc <- st_crop(cont, b)
  g <- g + geom_sf(data = zc, colour = "grey60", linewidth = 0.2) +
    geom_sf(data = st_crop(towers, b), shape = 24, fill = "black", colour = "white", size = 2.6) +
    geom_sf_text(data = st_crop(towers, b), aes(label = name), size = 2.4, nudge_y = -12) +
    geom_sf(data = t, aes(fill = species, shape = species), size = 2.4, colour = "grey15", stroke = 0.35) +
    scale_fill_manual(values = species_colors, name = NULL, drop = FALSE) +
    scale_shape_manual(values = spp_shapes, name = NULL, drop = FALSE) +
    coord_sf(xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"]), crs = st_crs(CRS_MA), datum = st_crs(4326), expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.6) +
    labs(tag = tag, title = title) + base_theme + theme(legend.position = "bottom", legend.text = element_text(size = 8, face = "italic"))
  g
}
pc <- zoom("Wetland (BGS)", "c", "Black Gum Swamp (wetland)")
pd <- zoom("Upland (EMS)", "d", "EMS tower footprint (upland)")

# ---- panel a: location ----
sdir <- file.path(tempdir(), "states"); unzip(file.path(S, "cb_2018_us_state_500k.zip"), exdir = sdir)
states <- st_read(list.files(sdir, pattern = "\\.shp$", full.names = TRUE)[1], quiet = TRUE) %>%
  filter(STUSPS %in% c("MA", "CT", "RI", "NH", "VT", "ME", "NY")) %>% st_transform(5070)
hf <- st_transform(st_as_sf(data.frame(lon = -72.18, lat = 42.537), coords = c("lon", "lat"), crs = 4326), 5070)
pa <- ggplot() + geom_sf(data = states, aes(fill = STUSPS == "MA"), colour = "grey40", linewidth = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#DCE6CC", `FALSE` = "grey95")) +
  geom_sf(data = hf, shape = 21, size = 3, fill = "#B03A2E", colour = "white", stroke = 0.6) +
  geom_sf_text(data = hf, label = "Harvard Forest", size = 2.6, nudge_y = 60000, fontface = "bold") +
  coord_sf(xlim = c(1.75e6, 2.15e6), ylim = c(2.25e6, 2.75e6), expand = FALSE, datum = NA) +
  labs(tag = "a") + base_theme + theme(axis.text = element_blank(), axis.ticks = element_blank())

p <- (pa | pb) / (pc | pd) + plot_layout(widths = c(1, 2.2), heights = c(1.25, 1))
ggsave(file.path(OUT, "site_map.png"), p, width = 10, height = 9.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "site_map.pdf"), p, width = 10, height = 9.5, bg = "white")
message("Saved site_map.png/pdf")
