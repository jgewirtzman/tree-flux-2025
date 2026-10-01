# ============================================================
# 08_site_map.R
# Study-site map (publication layout, 180 mm wide):
#   (a) Prospect Hill tract: shaded relief with elevation tint and 10-m contours, Black Gum
#       Swamp, the 2025 prism-survey grid, the study trees at both sites and the flux towers;
#       inset locates Harvard Forest in New England.
#   (b) Black Gum Swamp study trees by species.   (c) EMS study trees by species.
# Tree positions: field GPS surveys (BGS 2025-04-17, EMS 2026-01-16).
#
# Inputs: data/package/trees.csv (GPS positions), data/package/stand/ (BGS_VRP_2025.csv,
#   Black_Gum_Swamp.kmz); data/external/gis/ (elevation_ned/ 30-m DEM in MA State Plane,
#   elevation_contours.geojson, cb_2018_us_state_500k.zip; sources in data/README.md)
# Tower coordinates: AmeriFlux site metadata (US-Ha1, US-Ha2, US-xHA).
# Outputs: outputs/figures/site_map.png (600 dpi) / .pdf; outputs/tables/study_tree_coordinates.csv
# ============================================================
suppressPackageStartupMessages({
  library(sf); library(terra); library(ggplot2); library(dplyr); library(patchwork); library(ggspatial)
  library(ggnewscale); library(ggrepel)
})
S <- "data/package/stand"; GIS <- "data/external/gis"; OUT <- "outputs/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CRS_MA <- 26986
FONT <- "Helvetica"; BASE <- 8

species_colors <- c("N. sylvatica" = "#2A7F7A", "Q. rubra" = "#6E8B3D", "A. rubrum" = "#A7DAD1", "T. canadensis" = "#C9D6A4")
spp_shapes <- c("N. sylvatica" = 22, "A. rubrum" = 24, "T. canadensis" = 21, "Q. rubra" = 23)
site_colors <- c("Wetland (Black Gum Swamp)" = "#1F6F6A", "Upland (EMS)" = "#5E7A2E")
sp_name <- c(BG = "N. sylvatica", RM = "A. rubrum", HEM = "T. canadensis", RO = "Q. rubra")

# ---- data ----
trees <- read.csv(file.path("data", "package", "trees.csv")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  transmute(tree = Tree, species = factor(sp_name[toupper(species)], names(species_colors)),
            site = ifelse(plot == "BGS", "Wetland (Black Gum Swamp)", "Upland (EMS)")) %>%
  st_transform(CRS_MA)
stopifnot(nrow(trees) == 60, !anyNA(trees$species))
write.csv(cbind(st_drop_geometry(trees), st_coordinates(st_transform(trees, 4326))) %>% rename(lon = X, lat = Y),
          "outputs/tables/study_tree_coordinates.csv", row.names = FALSE)
towers <- st_as_sf(data.frame(name = c("EMS tower (US-Ha1)", "Hemlock tower (US-Ha2)", "NEON tower (US-xHA)"),
                              lon = c(-72.1715, -72.1779, -72.17265), lat = c(42.5378, 42.5393, 42.53691)),
                   coords = c("lon", "lat"), crs = 4326) %>% st_transform(CRS_MA)
kdir <- file.path(tempdir(), "bgs_kmz"); unzip(file.path(S, "Black_Gum_Swamp.kmz"), exdir = kdir)
swamp_line <- st_read(list.files(kdir, pattern = "\\.kml$", full.names = TRUE)[1], quiet = TRUE) %>% st_zm() %>%
  st_transform(CRS_MA) %>% st_union()
swamp <- st_collection_extract(st_polygonize(swamp_line), "POLYGON") %>% st_union()
vrp <- read.csv(file.path(S, "BGS_VRP_2025.csv")) %>% distinct(plot, lat, long) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% st_transform(CRS_MA)

# ---- panel a extent and relief ----
bb <- st_bbox(st_union(st_geometry(trees), st_geometry(towers)))
pad_x <- 380; pad_y <- 260
eb <- c(bb[["xmin"]] - pad_x, bb[["xmax"]] + pad_x, bb[["ymin"]] - pad_y, bb[["ymax"]] + pad_y)
dem <- rast(file.path(GIS, "elevation_ned")); crs(dem) <- paste0("EPSG:", CRS_MA)
dem_a <- crop(dem, ext(eb + c(-200, 200, -200, 200)))
dem_f <- disagg(dem_a, 3, method = "bilinear")                    # smoother relief at print size
sl <- terrain(dem_f, "slope", unit = "radians"); asp <- terrain(dem_f, "aspect", unit = "radians")
hs <- mean(rast(lapply(c(270, 315, 360), function(az) shade(sl, asp, angle = 40, direction = az))))   # multidirectional
hs_df <- as.data.frame(hs, xy = TRUE); names(hs_df)[3] <- "hs"
el_df <- as.data.frame(dem_f, xy = TRUE); names(el_df)[3] <- "elev"
cbb <- st_bbox(c(xmin = eb[1], xmax = eb[2], ymin = eb[3], ymax = eb[4]), crs = st_crs(CRS_MA))
cont <- st_read(file.path(GIS, "elevation_contours.geojson"), quiet = TRUE) %>% st_transform(CRS_MA) %>% st_crop(cbb)
cont$major <- (as.numeric(cont$ELEV) %% 50) == 0

theme_map <- theme_bw(base_size = BASE, base_family = FONT) +
  theme(axis.title = element_blank(), panel.grid = element_blank(), axis.text = element_text(size = BASE - 1.5, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.25), panel.border = element_rect(linewidth = 0.4, colour = "grey20"),
        plot.tag = element_text(face = "bold", size = BASE + 2), plot.title = element_text(size = BASE, face = "bold", hjust = 0),
        legend.text = element_text(size = BASE - 1), legend.title = element_text(size = BASE - 1),
        legend.key.size = unit(3.2, "mm"), legend.background = element_blank(), plot.margin = margin(2, 4, 2, 2))

tw_lab <- cbind(st_drop_geometry(towers), st_coordinates(towers))
sw_c <- st_coordinates(st_point_on_surface(swamp))
pa <- ggplot() +
  geom_raster(data = el_df, aes(x, y, fill = elev)) +
  scale_fill_gradientn(colours = c("#F4F2EA", "#E4E6D4", "#CCD5B6", "#ABBB94"), name = "Elevation (m)",
                       breaks = seq(320, 420, 20), guide = guide_colourbar(direction = "horizontal", barwidth = unit(24, "mm"), barheight = unit(2, "mm"), title.position = "top", order = 2)) +
  new_scale_fill() +
  geom_raster(data = hs_df, aes(x, y, alpha = hs), fill = "grey10", show.legend = FALSE) +
  scale_alpha_continuous(range = c(0.32, 0)) +
  geom_sf(data = cont %>% filter(!major), colour = "grey45", linewidth = 0.12, alpha = 0.55) +
  geom_sf(data = cont %>% filter(major), colour = "grey35", linewidth = 0.28, alpha = 0.7) +
  geom_sf(data = swamp, fill = "#7FB7B1", alpha = 0.35, colour = "#1F5E5A", linewidth = 0.45) +
  geom_sf(data = vrp, shape = 3, size = 0.55, stroke = 0.3, colour = "#1F5E5A", alpha = 0.8) +
  geom_sf(data = trees, aes(fill = site), shape = 21, colour = "white", stroke = 0.25, size = 1.6) +
  scale_fill_manual(values = site_colors, name = "Study trees", guide = guide_legend(order = 1, override.aes = list(size = 2))) +
  geom_sf(data = towers, shape = 24, fill = "black", colour = "white", size = 2.1, stroke = 0.4) +
  geom_text_repel(data = tw_lab, aes(X, Y, label = name), size = (BASE - 1.5) / .pt, family = FONT, box.padding = 0.5,
                  point.padding = 0.35, min.segment.length = 0, segment.size = 0.25, seed = 4, bg.color = "white", bg.r = 0.12) +
  annotate("text", x = sw_c[1], y = sw_c[2] - 170, label = "Black Gum\nSwamp", size = (BASE - 1) / .pt,
           family = FONT, fontface = "italic", colour = "#1F5E5A", lineheight = 0.9) +
  coord_sf(xlim = eb[1:2], ylim = eb[3:4], crs = st_crs(CRS_MA), datum = st_crs(4326), expand = FALSE) +
  scale_x_continuous(breaks = seq(-72.19, -72.16, by = 0.005)) + scale_y_continuous(breaks = seq(42.530, 42.545, by = 0.003)) +
  annotation_scale(location = "bl", width_hint = 0.18, text_cex = 0.6, height = unit(1.4, "mm"), text_family = FONT, line_width = 0.4) +
  annotation_north_arrow(location = "tr", height = unit(6, "mm"), width = unit(4.5, "mm"), style = north_arrow_orienteering(text_size = 5)) +
  labs(tag = "a") + theme_map +
  theme(legend.position = "inside", legend.position.inside = c(0.995, 0.04), legend.justification = c(1, 0),
        legend.box = "horizontal", legend.box.just = "bottom",
        legend.box.background = element_rect(fill = alpha("white", 0.85), colour = NA), legend.margin = margin(2, 3, 2, 3))

# ---- inset locator ----
sdir <- file.path(tempdir(), "states"); unzip(file.path(GIS, "cb_2018_us_state_500k.zip"), exdir = sdir)
states <- st_read(list.files(sdir, pattern = "\\.shp$", full.names = TRUE)[1], quiet = TRUE) %>%
  filter(STUSPS %in% c("MA", "CT", "RI", "NH", "VT", "ME", "NY")) %>% st_transform(5070)
hf <- st_transform(st_as_sf(data.frame(lon = -72.18, lat = 42.537), coords = c("lon", "lat"), crs = 4326), 5070)
inset <- ggplot() + geom_sf(data = states, aes(fill = STUSPS == "MA"), colour = "grey45", linewidth = 0.15, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#DCE6CC", `FALSE` = "white")) +
  geom_sf(data = hf, shape = 21, size = 1.6, fill = "#B03A2E", colour = "white", stroke = 0.3) +
  coord_sf(xlim = c(1.78e6, 2.12e6), ylim = c(2.30e6, 2.70e6), expand = FALSE, datum = NA) +
  theme_void() + theme(panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.3),
                       panel.background = element_rect(fill = "white", colour = NA))
pa <- pa + inset_element(inset, left = 0.006, bottom = 0.6, right = 0.185, top = 0.994, align_to = "panel")

# ---- zoom panels ----
zoom <- function(site_name, tag, title, pad = 18, legend = FALSE) {
  t <- trees %>% filter(site == site_name)
  b <- st_bbox(st_buffer(st_union(t), pad))
  g <- ggplot()
  if (grepl("Wetland", site_name)) g <- g + geom_sf(data = st_crop(swamp, b), fill = "#7FB7B1", alpha = 0.25, colour = "#1F5E5A", linewidth = 0.4)
  g + geom_sf(data = st_crop(cont, b), colour = "grey55", linewidth = 0.18) +
    geom_sf(data = st_crop(towers, b), shape = 24, fill = "black", colour = "white", size = 2.1, stroke = 0.4) +
    geom_sf(data = t, aes(fill = species, shape = species), size = 2, colour = "grey15", stroke = 0.3) +
    scale_fill_manual(values = species_colors, name = NULL, limits = names(species_colors), drop = FALSE) +
    scale_shape_manual(values = spp_shapes, name = NULL, limits = names(species_colors), drop = FALSE, guide = "none") +
    guides(fill = guide_legend(override.aes = list(shape = unname(spp_shapes[names(species_colors)]),
                                                   fill = unname(species_colors), size = 2.2, colour = "grey15"))) +
    annotate("text", x = -Inf, y = Inf, label = title, hjust = -0.06, vjust = 1.6, size = BASE / .pt, family = FONT, fontface = "bold") +
    coord_sf(xlim = c(b[["xmin"]], b[["xmax"]]), ylim = c(b[["ymin"]], b[["ymax"]]), crs = st_crs(CRS_MA), datum = NA, expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.6, height = unit(1.2, "mm"), text_family = FONT, line_width = 0.4) +
    labs(tag = tag) + theme_map +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.text = element_text(size = BASE, face = "italic"),
          legend.position = if (legend) "bottom" else "none")
}
pb <- zoom("Wetland (Black Gum Swamp)", "b", "Black Gum Swamp")
neon <- cbind(st_drop_geometry(towers), st_coordinates(towers)) %>% filter(grepl("NEON", name))
e5 <- trees %>% filter(site == "Upland (EMS)") %>% mutate(x = st_coordinates(.)[, 1]) %>% filter(x < min(x) + 40)
e5c <- st_coordinates(st_centroid(st_union(e5)))
pc <- zoom("Upland (EMS)", "c", "EMS tower footprint") +
  annotate("text", x = e5c[1] + 14, y = e5c[2] - 14, label = "plot E5", size = (BASE - 1.5) / .pt, family = FONT, hjust = 0, colour = "grey30") +
  geom_text_repel(data = neon, aes(X, Y, label = "NEON tower"), size = (BASE - 1.5) / .pt, family = FONT, nudge_y = -22, segment.size = 0.2, seed = 2)

# one species legend for panels b and c, built from all four species
legd <- data.frame(x = 1:4, y = 1, species = factor(names(species_colors), names(species_colors)))
leg <- cowplot::get_plot_component(ggplot(legd, aes(x, y, fill = species, shape = species)) +
  geom_point(size = 2.2, colour = "grey15", stroke = 0.3) +
  scale_fill_manual(values = species_colors, name = NULL) + scale_shape_manual(values = spp_shapes, name = NULL) +
  theme_void(base_size = BASE, base_family = FONT) +
  theme(legend.position = "bottom", legend.text = element_text(size = BASE, face = "italic"), legend.key.size = unit(3.5, "mm")),
  "guide-box-bottom")
p <- pa / (pb + pc + plot_layout(widths = c(1, 2.2))) / patchwork::wrap_elements(leg) + plot_layout(heights = c(1.55, 1, 0.08))
ggsave(file.path(OUT, "site_map.pdf"), p, width = 180, height = 165, units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(OUT, "site_map.png"), p, width = 180, height = 165, units = "mm", dpi = 600, bg = "white")
message("Saved site_map.png/pdf")

