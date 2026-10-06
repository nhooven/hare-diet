# PROJECT: Diet
# SCRIPT: 14 - REVIEW - Browse ranks
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 06 Oct 2026
# COMPLETED: 
# LAST MODIFIED: 
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# from studies
data.browse <- read.csv("data_review/data_browse.csv")

# ______________________________________________________________________________
# 3. Preliminary plot ----

# let's just see if any patterns jump out

# ______________________________________________________________________________

ggplot(data = data.browse) +
  
  theme_classic() +
  
  geom_point(aes(x = importance_rank,
                 y = reorder(Genus, 1 / importance_rank)),
             shape = 21) +
  
  scale_x_reverse() +
  
  theme(axis.title = element_blank())

# and just the taxa we have
loomis.genera <- c("Pinus", "Larix", "Ribes", "Rosa", "Pseudotsuga",
                   "Vaccinium", "Salix", "Picea", "Abies", "Amelanchier",
                   "Populus", "Alnus", "Arctostaphylos", "Lonicera",
                   "Paxistima")

ggplot(data = data.browse |> filter(Genus %in% loomis.genera)) +
  
  theme_classic() +
  
  geom_point(aes(x = importance_rank,
                 y = reorder(Genus, 1 / importance_rank)),
             shape = 21) +
  
  scale_x_reverse() +
  
  theme(axis.title = element_blank())

# NEXT:
  # a way to re-calculate ranks after filtering?
  # a mean rank calculation