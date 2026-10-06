# PROJECT: Diet
# SCRIPT: 12 - REVIEW - Map
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 06 Oct 2026
# COMPLETED: 
# LAST MODIFIED: 06 Oct 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)
library(geobounds)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# studies and coordinates
data.map <- read.csv("data_review/data_map.csv") |>
  
  st_as_sf(coords = c("long", "lat"),
           crs = "epsg:4326")

# hare range map
all.ranges <- st_read("data_review/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")

hare.range <- all.ranges |> filter(sci_name == "Lepus americanus")

# remove
rm(all.ranges)

plot(st_geometry(hare.range))

# write to file
st_write(hare.range, "data_review/hare_range.shp")

# North American political boundaries
na.bounds <- gb_get_adm1(country = c("United States", "Canada"),
                         simplified = T)

plot(st_geometry(na.bounds))

# write to file
st_write(na.bounds, "data_review/na_bounds.shp", append = F)

# ______________________________________________________________________________
# 3. Plot ----
# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  # NA political boundaries
  geom_sf(data = na.bounds,
          fill = "white") +
  
  # hare IUCN range
  geom_sf(data = hare.range,
          fill = "dodgerblue3",
          alpha = 0.5) +
  
  # hare studies
  geom_sf(data = data.map,
          aes(shape = aspect),
          fill = "white") +
  
  # map frame
  coord_sf(xlim = c(-170, -52),
           ylim = c(35, 70),
           expand = F) +
  
  # theme
  theme(
    
    # legend
    legend.position = c(0.2, 0.2),
    legend.title = element_blank(),
    legend.background = element_rect(color = "black"),
    
    # axis ticks
    axis.text = element_blank(),
    axis.ticks = element_blank()
    
  ) +
  
  # shapes
  scale_shape_manual(values = c(21, 22, 23))



# NEXT:
  # manually jitter the Kluane studies
  # change size/shapes (colors?)
  # come up with a better word than "aspect"
  # fix legend spacing