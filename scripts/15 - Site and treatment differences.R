# PROJECT: Diet
# SCRIPT: 15 - Site and treatment differences
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 07 Oct 2026
# COMPLETED: 07 Oct 2026
# LAST MODIFIED: 07 Oct 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(vegan)    # NMDS

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# taxa information (from script 07)
reads.taxa <- read.csv("data_cleaned/reads_cleaned.csv")

# ______________________________________________________________________________
# 3. Clean ----

# vegan accepts a n x variable "community data matrix"
# so samples need to be rows, and taxa are columns

# I want a presence/absence and a read count matrix

# ______________________________________________________________________________
# 3a. Function to prep ----
# ______________________________________________________________________________

prep_cdm <- function (.reads,
                      type = c("po", "read")) {
  
  if (type == "po") {
    
    # data.frame
    taxa.df <- .reads |>
      
      # keep only necessary columns (sample, Final.taxon)
      dplyr::select(Final.sample.ID, final.taxon) |>
      
      # add a presence column
      mutate(presence = 1) |>
      
      # cast wide
      pivot_wider(names_from = Final.sample.ID,
                  values_from = presence) |>
      
      # coerce everthing to a character
      mutate(across(everything(), as.character)) |>
      
      # replace "NULL" with 0, all else 1
      mutate(across(starts_with("S"), ~ifelse(. == "NULL",
                                              0,
                                              1))) |>
      
      # coerce to Boolean
      mutate(across(starts_with("S"), as.integer))
    
    # remove first column and transpose
    cdm <- taxa.df[ , -1] |> as.matrix() |> t()
    
    # else, use reads
  } else {
    
    # data.frame
    taxa.df <- .reads |>
      
      # keep only necessary columns (sample, Final.taxon)
      dplyr::select(Final.sample.ID, final.taxon, reads) |>
      
      # sum reads for the sample taxon
      group_by(Final.sample.ID, final.taxon) |>
      summarize(reads = sum(reads)) |>
      
      # cast wide
      pivot_wider(names_from = Final.sample.ID,
                  values_from = reads)
    
    # remove first column and transpose and replace NAs
    cdm <- taxa.df[ , -1] |> as.matrix() |> t()
    
    cdm[which(is.na(cdm))] <- 0
    
  }
  
  return(cdm)
  
}

# ______________________________________________________________________________
# 3b. Use function ----

reads.off <- reads.taxa |> filter(Season == "Off")
reads.on <- reads.taxa |> filter(Season == "On")

# ______________________________________________________________________________

cdm.off.po <- prep_cdm(reads.off, "po")
cdm.off.read <- prep_cdm(reads.off, "read")

cdm.on.po <- prep_cdm(reads.on, "po")
cdm.on.read <- prep_cdm(reads.on, "read")

# ______________________________________________________________________________
# 4. Fit NMDS ----
# ______________________________________________________________________________
# 4a. Wrapper function ----
# ______________________________________________________________________________

fit_nmds <- function (.cdm) {
  
  nmds.out <- metaMDS(
    
    comm = .cdm,
    distance = "bray",
    k = 2,
    autotransform = F,
    try = 40
    
  )
  
  return(nmds.out)
  
}

# ______________________________________________________________________________
# 4b. Fit ----
# ______________________________________________________________________________

nmds.off.po <- fit_nmds(cdm.off.po)
nmds.off.read <- fit_nmds(cdm.off.read)
nmds.on.po <- fit_nmds(cdm.on.po)
nmds.on.read <- fit_nmds(cdm.on.read)

# ______________________________________________________________________________
# 5. Evaluate ----
# ______________________________________________________________________________

# stress
nmds.off.po$stress
nmds.off.read$stress
nmds.on.po$stress
nmds.on.read$stress

# GOF
hist(goodness(nmds.off.po))
hist(goodness(nmds.off.read))
hist(goodness(nmds.on.po))
hist(goodness(nmds.on.read))

# distance vs. dissimilarity
plot(nmds.off.po$diss, nmds.off.po$dist)
plot(nmds.off.read$diss, nmds.off.read$dist)
plot(nmds.on.po$diss, nmds.on.po$dist)
plot(nmds.on.read$diss, nmds.on.read$dist)

# ______________________________________________________________________________
# 6. Extract scores for plotting ----
# ______________________________________________________________________________
# 6a. Function ----
# ______________________________________________________________________________

extract_nmds <- function (.reads,
                          .nmds) {
  
  extract.out <- .reads |>
    
    dplyr::select(Final.sample.ID,
                  Treatment,
                  Season,
                  Sex,
                  Site) |>
    
    # keep only one row per sample
    group_by(Final.sample.ID) |>
    slice(1) |>
    ungroup() |>
    
    # join in site-specific NMDS scores
    left_join(
      
      data.frame(
        
        Final.sample.ID = rownames(scores(.nmds)$sites),
        NMDS1 = scores(.nmds)$sites[ , 1],
        NMDS2 = scores(.nmds)$sites[ , 2]
        
      )
      
    )
  
  return(extract.out)
  
}

# ______________________________________________________________________________
# 6b. Use, add in identifiers for plotting ----
# ______________________________________________________________________________

all.nmds.scores <- bind_rows(
  
  extract_nmds(reads.off, nmds.off.po) |> mutate(metric = "po"),
  extract_nmds(reads.off, nmds.off.read) |> mutate(metric = "read"),
  extract_nmds(reads.on, nmds.on.po) |> mutate(metric = "po"),
  extract_nmds(reads.on, nmds.on.read) |> mutate(metric = "read")
  
) |>
  
  # site cluster
  mutate(clust = substr(Site, 1, 1)) |>
  
  # factor levels/labels
  mutate(
    
    Treatment = factor(Treatment, levels = c("Control", "Retention", "Piling")),
    Season = factor(Season, labels = c("snow-off", "snow-on")),
    metric = factor(metric, labels = c("presence/absence", "reads"))
    
  )

# ______________________________________________________________________________
# 7. Plots ----
# ______________________________________________________________________________
# 7a. Sex ----
# ______________________________________________________________________________

ggplot(data = all.nmds.scores) +
  
  theme_bw() +
  
  facet_grid(metric ~ Season) +
  
  # ellipses
  stat_ellipse(aes(x = NMDS1,
                   y = NMDS2,
                   linetype = Sex,
                   color = Sex)) +
  
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 color = Sex,
                 shape = Sex),
             size = 1.1) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA,
                                        fill = "lightgray"),
        strip.text = element_text(hjust = 0),
        legend.title = element_blank()) +
  
  # colors, shapes, etc.
  scale_color_manual(values = c("#FF3300", "gray50"))

# 495 x 414

# ______________________________________________________________________________
# 7b. Treatment ----
# ______________________________________________________________________________

ggplot(data = all.nmds.scores) +
  
  theme_bw() +
  
  facet_grid(metric ~ Season) +
  
  # ellipses
  stat_ellipse(aes(x = NMDS1,
                   y = NMDS2,
                   linetype = Treatment,
                   color = Treatment)) +
  
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 color = Treatment,
                 shape = Treatment),
             size = 1.1) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA,
                                        fill = "lightgray"),
        strip.text = element_text(hjust = 0),
        legend.title = element_blank()) +
  
  # colors, shapes, etc.
  scale_color_manual(values = c("gray", "purple", "orange"))

# 495 x 414

# ______________________________________________________________________________
# 7c. Cluster ----
# ______________________________________________________________________________

ggplot(data = all.nmds.scores) +
  
  theme_bw() +
  
  facet_grid(metric ~ Season) +
  
  # ellipses
  stat_ellipse(aes(x = NMDS1,
                   y = NMDS2,
                   linetype = clust,
                   color = clust)) +
  
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 color = clust,
                 shape = clust),
             size = 1.1) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA,
                                        fill = "lightgray"),
        strip.text = element_text(hjust = 0),
        legend.title = element_blank()) +
  
  # colors, shapes, etc.
  scale_color_manual(values = c("darkgreen", "green3", "gold"))

# 495 x 414


# NEXT:
  # PERMANOVA
  # ensure that treatment assignment is correct?