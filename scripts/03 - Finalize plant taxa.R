# PROJECT: Diet
# SCRIPT: 03 - Finalize plant taxa
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 08 Sep 2026
# COMPLETED: 08 Sep 2026
# LAST MODIFIED: 08 Sep 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 0. Purpose ----
# ______________________________________________________________________________

# after subsequent BLASTs, we should group all relevant identified ASVs into one
# table, so we can attribute them to sample reads in the next script

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

asvs.1 <- read.csv("data_cleaned/asvs_identified.csv")
asvs.2 <- read.csv("data_cleaned/asvs_identified2.csv")
asvs.3 <- read.csv("data_cleaned/asvs_identified3.csv")
asvs.4 <- read.csv("data_cleaned/asvs_identified4.csv")

# ______________________________________________________________________________
# 3. Clean and bind together ----
# ______________________________________________________________________________

# get original ID file into the right format
str(asvs.1)
str(asvs.2)

asvs.1 <- asvs.1 |>
  
  dplyr::select(ASV.unq, batchID, ASV, 
                n_reads, cover, ident, 
                taxon_inferred, family, genus, species) |>
  
  rename(taxa = taxon_inferred)

# keep columns for other ASV files
asvs.2 <- asvs.2 |> dplyr::select(ASV.unq, batchID, ASV,
                                  n_reads, cover, ident,
                                  taxa, family, genus, species)

asvs.3 <- asvs.3 |> dplyr::select(ASV.unq, batchID, ASV,
                                  n_reads, cover, ident,
                                  taxa, family, genus, species)

asvs.4 <- asvs.4 |> dplyr::select(ASV.unq, batchID, ASV,
                                  n_reads, cover, ident,
                                  taxa, family, genus, species)

# bind in
asvs.all <- bind_rows(asvs.1, asvs.2, asvs.3, asvs.4)

# remove contaminant reads
asvs.all <- asvs.all[-which(asvs.all$genus == "Glycine"), ]

# ______________________________________________________________________________
# 4. Add missing taxonomic information ----
# ______________________________________________________________________________

# missing family
(miss.fam <- which(asvs.all$family == ""))

asvs.all$genus[miss.fam]

asvs.all$family[miss.fam] <- c("Ranunculaceae", "Grossulariaceae", 
                               "Salicaceae", "Cupressaceae",
                               "Adoxaceae")

# missing species
# Aquilegia flavescens
asvs.all$genus[7]
asvs.all$species[7] <- "flavescens"

# Galium triflorum
asvs.all$genus[187]
asvs.all$genus[227]
asvs.all$species[187] <- "triflorum"
asvs.all$species[227] <- "triflorum"

# Viola adunca
asvs.all$genus[274]
asvs.all$species[274] <- "adunca"

# ______________________________________________________________________________
# 5. Final taxon name ----
# ______________________________________________________________________________

asvs.all.1 <- asvs.all |>
  
  mutate(
    
    final.taxon = case_when(
      
      # just a family
      genus == "" & species == "" ~ family,
      
      # just a genus
      genus != "" & species == "" ~ genus,
      
      # species included
      genus != "" & species != "" ~ paste0(genus, " ", species)
      
    )
    
  ) |>
  
  # drop "taxa" notes column
  dplyr::select(-taxa)

# ______________________________________________________________________________
# 6. How many unique taxa? ----
# ______________________________________________________________________________

length(unique(asvs.all.1$final.taxon))

# ______________________________________________________________________________
# 7. Write to file ----
# ______________________________________________________________________________

write.csv(asvs.all.1, "data_cleaned/all_taxa.csv", row.names = F)
