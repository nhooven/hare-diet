# PROJECT: Diet
# SCRIPT: 18 - Study site map
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 16 Oct 2026
# COMPLETED: 16 Oct 2026
# LAST MODIFIED: 16 Oct 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(cowplot)
library(ggspatial)       # scale bar

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# paths
dir.spat <- "D:/hare_project/data_spatial/"
dir.ch2 <- "D:/hare_project/Maps/Ch 2 study area map/rasters/"

# WA boundary
wa.border <- st_read(paste0(dir.spat, "Borders/better_WA_border/wa_border_1.shp"))

# Loomis boundary
lsf <- st_read(paste0(dir.spat, "Loomis layers/LSF boundary/loomis_boundary_wgs.shp"))

# PCT units
pct.units <- st_read(paste0(dir.spat, "Units/units_fixed_utm/units_fixed_utm.shp"))

# hillshade
hillshade <- rast(paste0(dir.ch2, "hs_lsf.tif"))

# canopy cover
cc <- rast(paste0(dir.spat, "Rasters/TCC/tcc_2024_utm.tif"))

# ______________________________________________________________________________
# 3. Process spatial layers ----
# ______________________________________________________________________________

# UTM crs
crs.utm <- "epsg:32611"

# aggregate hillshade
hillshade.1 <- terra::aggregate(hillshade, fact = 30)

# crop and mask CC
cc.1 <- cc |> project(crs.utm) |>
  
  crop(st_transform(lsf, crs(cc))) |>
  
  mask(st_transform(lsf, crs(cc))) |>
  
  # clamp to 100
  clamp(upper = 100)

# PCT units
pct.units <- pct.units |> dplyr::select(name, geometry) |>
  
  filter(name %in% c("1A", "1B", "1C",
                     "2A", "2B", "2C",
                     "3A", "3B", "3C")) |>
  
  arrange(name)

# cluster centroids
unit.cent <- bind_rows(
  
  pct.units |> filter(name %in% c("1A", "1B", "1C")) |> st_combine() |> st_centroid() |> st_as_sf(),
  pct.units |> filter(name %in% c("2A", "2B", "2C")) |> st_combine() |> st_centroid() |> st_as_sf(),
  pct.units |> filter(name %in% c("3A", "3B", "3C")) |> st_combine() |> st_centroid() |> st_as_sf()
  
)

# ______________________________________________________________________________
# 4. Map with insets ----
# ______________________________________________________________________________
# 4a. Washington with LSF boundary ----

# washington centroid
wa.cent <- st_coordinates(st_centroid(wa.border))

# ______________________________________________________________________________

wa.map <- ggplot() +
  
  theme_bw() +
  
  geom_sf(data = wa.border,
          color = NA,
          fill = "gray85") +
  
  # LSF boundary
  geom_sf(data = lsf,
          color = NA,
          fill = "darkgreen") +
  
  coord_sf(datum = st_crs(32611)) +
  
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank()) +
  
  annotate(geom = "text",
           x = wa.cent[1],
           y = wa.cent[2],
           label = "WA")

wa.map

# ______________________________________________________________________________
# 4b. Full LSF ----
# ______________________________________________________________________________

lsf.map <- ggplot() +
  
  theme_bw() +
  
  # cc raster
  geom_raster(data = cc.1,
              aes(x = x,
                  y = y,
                  fill = tcc_2024_utm)) +
  
  # fill
  scale_fill_gradient(low = "white",
                      high = "darkgreen",
                      na.value = "white") +
  
  # unit centroids
  geom_sf(data = unit.cent,
          fill = "gold",
          stroke = 0.75,
          size = 2,
          shape = 23) + 
  
  # LSF boundary
  geom_sf(data = lsf,
          fill = NA) +
  
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank()) +
  
  # scale bar
  annotation_scale(aes(location = "bl",
                       style = "ticks"),
                   line_width = 0.5,
                   height = unit(0.1, "cm"),
                   text_cex = 0.6,
                   pad_x = unit(0.2, "cm"),
                   pad_y = unit(0.2, "cm"),
                   line_col = "black",
                   text_col = "black") +
  
  coord_sf(datum = st_crs(32611),
           expand = F)

lsf.map

# ______________________________________________________________________________
# 4c. Plot all ----
# ______________________________________________________________________________

ggdraw(lsf.map) +
  
  draw_plot(wa.map,
            x = 0.23,
            y = 0.62,
            width = 0.2,
            height = 0.2)

ggsave(filename = "other_plots/LSF_map.png",
       dpi = 300,
       width = 574,
       height = 380,
       scale = 3.5,
       units = "px",
       bg = "white")
