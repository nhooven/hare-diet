# PROJECT: Diet
# SCRIPT: 03 - Finalize plant taxa
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 08 Sep 2026
# COMPLETED: 08 Sep 2026
# LAST MODIFIED: 14 Sep 2026
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

# Asteraceae sub-groups
# add column
asvs.all$sub.group <- ""

(unk.Aster <- which(asvs.all$family == "Asteraceae" & asvs.all$genus == ""))

asvs.all$taxa[unk.Aster]

# change genus
asvs.all$genus[unk.Aster][1:4] <- c("Solidago / Symphyotrichum", "Hieracium / Taraxacum",
                                    "Solidago / Symphyotrichum", "Hieracium / Taraxacum")

# assign sub-group
asvs.all$sub.group[unk.Aster] <- c("Asteraceae 1", "Asteraceae 2",
                                   "Asteraceae 1", "Asteraceae 2",
                                   "unknown Asteraceae")

# ______________________________________________________________________________
# 5. Final taxon name and functional groups ----
# ______________________________________________________________________________

asvs.all.1 <- asvs.all |>
  
  mutate(
    
    final.taxon = case_when(
      
      # just a family (no sub-group)
      genus == "" & species == "" & sub.group == "" ~ family,
      
      # just a genus (no sub-group)
      genus != "" & species == "" & sub.group == "" ~ genus,
      
      # species included
      genus != "" & species != "" ~ paste0(genus, " ", species),
      
      # just a family and sub-group
      family != "" & sub.group != "" ~ sub.group
      
    )
    
  ) |>
  
  # assign functional groups
  # cat1 - broadest category
    # conifer
    # woody broadleaf
    # sub-shrub
    # forb
    # graminoid
    # unknown
  mutate(
    
    cat1 = case_when(
      
      genus %in% c("Pinus",
                   "Larix",
                   "Picea",
                   "Abies",
                   "Juniperus",
                   "Pseudotsuga") ~ "conifer",
      genus %in% c("Shepherdia",
                   "Alnus",
                   "Rubus",
                   "Ribes",
                   "Populus",
                   "Artemisia",
                   "Salix",
                   "Lonicera",
                   "Philadelphus",
                   "Rhododendron",
                   "Rosa",
                   "Celtis",
                   "Symphoricarpos") ~ "woody broadleaf",
      
      genus %in% c("Vaccinium",
                   "Paxistima",
                   "Arctostaphylos",
                   "Spiraea",
                   "Linnaea") ~ "sub-shrub",
      
      genus %in% c("Aquilegia",
                   "Lupinus",
                   "Trifolium",
                   "Chamaenerion",
                   "Conrus",
                   "Astragalus",
                   "Epilobium",
                   "Veronica",
                   "Viola",
                   "Thalictrum",
                   "Antennaria",
                   "Arnica",
                   "Geum",
                   "Verbena",
                   "Salvia",
                   "Eriogonum",
                   "Allium",
                   "Goodyera",
                   "Brassica",
                   "Rumex",
                   "Lathyrus",
                   "Ligusticum",
                   "Penstemon",
                   "Spergularia",
                   "Polygonum",
                   "Hedysarum",
                   "Ranunculus",
                   "Oenothera",
                   "Sedum",
                   "Amaranthus",
                   "Orthilia",
                   "Verbascum",
                   "Streptopus",
                   "Solidago / Symphyotrichum",
                   "Achillea",
                   "Adenocaulon",
                   "Convolvulus",
                   "Celtis",
                   "Lycium",
                   "Hieracium / Taraxacum",
                   "Erigeron",
                   "Anaphalis",
                   "Athyrium",
                   "Descurainia",
                   "Chorispora",
                   "Lactuca",
                   "Capsella",
                   "Arenaria",
                   "Chenopodium",
                   "Gnaphalium",
                   "Osmorhiza",
                   "Cerastium",
                   "Cornus",
                   "Galium",
                   "Stellaria") ~ "forb",
      
      family %in% c("Poaceae", "Cyperaceae") ~ "graminoid",
      
      genus %in% c("") ~ "unknown"
      
    )
    
  ) |>
  
  # cat2 - adds a lodgepole level
  mutate(
    
    cat2 = ifelse(species == "contorta",
                  "lodgepole pine",
                  cat1)
    
  )

# add more as needed

# ______________________________________________________________________________
# 6. How many unique taxa? ----
# ______________________________________________________________________________

length(unique(asvs.all.1$final.taxon))

# by functional group
asvs.all.1 |> group_by(final.taxon) |> 
  
  slice(1) |>
  ungroup() |>
  group_by(cat1, cat2) |>
  summarize(n())

# ______________________________________________________________________________
# 7. Write to file ----
# ______________________________________________________________________________

write.csv(asvs.all.1, "data_cleaned/all_taxa.csv", row.names = F)
