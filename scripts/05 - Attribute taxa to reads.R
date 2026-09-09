# PROJECT: Diet
# SCRIPT: 05 - Attribute taxa to reads
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 09 Sep 2026
# COMPLETED: 09 Sep 2026
# LAST MODIFIED: 09 Sep 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

samples.lookup <- read.csv("data_cleaned/samples_lookup.csv")
reads.all <- read.csv("data_cleaned/reads_all.csv")
taxa <- read.csv("data_cleaned/all_taxa.csv")

# ______________________________________________________________________________
# 3. Join ----
# ______________________________________________________________________________

reads.taxa <- reads.all |>
  
  left_join(
    
    taxa |>
      
      dplyr::select(ASV.unq,
                    family,
                    genus,
                    species,
                    final.taxon),
    
    by = "ASV.unq"
    
  ) |>
  
  # remove NAs (bait, contaminants, unidentified ASVs)
  drop_na(final.taxon)

# ______________________________________________________________________________
# 4. Write ----
# ______________________________________________________________________________

write.csv(reads.taxa, "data_cleaned/reads_taxa.csv", row.names = F)
