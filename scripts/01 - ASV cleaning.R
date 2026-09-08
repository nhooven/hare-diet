# PROJECT: Diet
# SCRIPT: 01 - ASV cleaning
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 04 Sep 2026
# COMPLETED: 04 Sep 2026
# LAST MODIFIED: 04 Sep 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 0. Purpose ----
# ______________________________________________________________________________

# Here we'll:
  # read in the full ASV table (from three metabarcoding "batches")
  # sort by identification status from initial bioinformatics
  # remove obvious contaminant ASVs
  # extract ASVs to BLAST again
  # write FASTA files

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(seqinr)

# ______________________________________________________________________________
# 2. Read data ----
# ______________________________________________________________________________

asvs <- read.csv("data_raw/asvs.csv")

# ______________________________________________________________________________
# 3. Raw summaries ----
# ______________________________________________________________________________

# number of ASVs in each group
asv.group.summary <- asvs |>
  
  group_by(group) |>
  
  summarize(n.asvs = n(),
            n.reads = sum(n_reads))

asv.group.summary

# ______________________________________________________________________________
# 4. Cleaning ----

# group identifier
unique(asvs$group)

# ______________________________________________________________________________

asvs.1 <- asvs |>
  
  # unique ASV column
  mutate(ASV.unq = paste0(ASV, "_", batchID)) |>
  
  # remove contaminant and control ASVs
  filter(group %in% c("identified",
                      "grass",
                      ""))

# split
asvs.identified <- asvs.1 |> filter(group == "identified")
asvs.grass <- asvs.1 |> filter(group == "grass")
asvs.unidentified <- asvs.1 |> filter(group == "")

# ______________________________________________________________________________
# 5. Save to file ----
# ______________________________________________________________________________

write.csv(asvs.identified, "data_cleaned/asvs_identified.csv", row.names = F)
write.csv(asvs.grass, "data_cleaned/asvs_grass.csv", row.names = F)
write.csv(asvs.unidentified, "data_cleaned/asvs_unidentified.csv", row.names = F)

# ______________________________________________________________________________
# 6. Write FASTA files ----
# ______________________________________________________________________________

# grass
write.fasta(
  
  as.list(asvs.grass$sequence),
  asvs.grass$ASV.unq,
  "data_cleaned/FASTA/grass.fasta",
  as.string = T
  
)

# unidentified
write.fasta(
  
  as.list(asvs.unidentified$sequence),
  asvs.unidentified$ASV.unq,
  "data_cleaned/FASTA/unidentified.fasta",
  as.string = T
  
)
