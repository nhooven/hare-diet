# PROJECT: Diet
# SCRIPT: 12 - REVIEW - Map
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 06 Oct 2026
# COMPLETED: 12 Oct 2026
# LAST MODIFIED: 12 Oct 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)
library(geobounds)

# ______________________________________________________________________________
# 2. Spatial data ----
# ______________________________________________________________________________

# hare range map
#all.ranges <- st_read("data_review/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")

#hare.range <- all.ranges |> filter(sci_name == "Lepus americanus")

# remove
#rm(all.ranges)

#plot(st_geometry(hare.range))

# write to file
#st_write(hare.range, "data_review/hare_range.shp")

# North American political boundaries
#na.bounds <- gb_get_adm1(country = c("United States", "Canada"),
#                         simplified = T)

#plot(st_geometry(na.bounds))

# write to file
#st_write(na.bounds, "data_review/na_bounds.shp", append = F)

# read
hare.range <- st_read("data_review/hare_range.shp")
na.bounds <- st_read("data_review/na_bounds.shp")

# ______________________________________________________________________________
# 3. Studies and coordinates ----
# ______________________________________________________________________________

data.map <- read.csv("data_review/data_map.csv") |>
  
  st_as_sf(coords = c("long", "lat"),
           crs = "epsg:4326") |>
  
  # 
  
  # bind in our study for funs
  bind_rows(
    
    data.frame(
      
      Study = "Hooven et al.",
      Location = "Loomis, WA",
      
      # Crazy Beetle parking spot
      long = -119.813253,
      lat = 48.865928,
      
      loc = "exact",
      season = "both",
      aspect = "this study"
      
    ) |>
    
    st_as_sf(coords = c("long", "lat"),
             crs = "epsg:4326")
    
  ) |>
  
  # factor levels
  mutate(aspect = factor(aspect, levels = c("richness",
                                            "browse",
                                            "both",
                                            "this study")))

# jitter Kluane points
set.seed(415)

data.map.kluane <- data.map |> filter(Location == "Kluane Lake, YK") |> st_jitter(amount = 2)

data.map.rest <- data.map |> filter(Location != "Kluane Lake, YK")

# ______________________________________________________________________________
# 4. Plot ----
# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  # NA political boundaries
  geom_sf(data = na.bounds,
          color = "black",
          fill = "gray95",
          linewidth = 0.2) +
  
  # hare IUCN range
  geom_sf(data = hare.range,
          fill = "dodgerblue3",
          color = NA,
          alpha = 0.5) +
  
  # hare studies - Kluane
  geom_sf(data = data.map.kluane,
          aes(shape = aspect,
              size = aspect,
              fill = aspect,
              color = aspect,
              stroke = aspect)) +

  # hare studies 
  geom_sf(data = data.map.rest,
          aes(shape = aspect,
              size = aspect,
              fill = aspect,
              color = aspect,
              stroke = aspect)) +
  
  # map frame
  coord_sf(xlim = c(-170, -52),
           ylim = c(35, 70),
           expand = F) +
  
  # theme
  theme(
    
    # legend
    legend.position = c(0.2, 0.2),
    legend.title = element_blank(),
    legend.text = element_text(size = 7,
                               vjust = 0.5),
    legend.background = element_rect(color = "black"),
    legend.key.spacing = unit(0.1, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width = unit(0.35, "cm"),
    
    # axis ticks
    axis.text = element_blank(),
    axis.ticks = element_blank()
    
  ) +
  
  # shapes, sizes, colors
  scale_shape_manual(values = c(22, 23, 21, 4)) +
  scale_size_manual(values = c(2.1, 2, 2, 1.7)) +
  scale_fill_manual(values = c("green4", "dodgerblue3", "darkorange", "white")) +
  scale_color_manual(values = c(rep("white", 3), "black")) +
  scale_discrete_manual(values = c(0.5, 0.5, 0.5, 0.7),
                        aesthetics = "stroke")

# 574 x 380
