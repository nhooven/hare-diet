# PROJECT: Diet
# SCRIPT: 02 - Process BLAST results
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 08 Sep 2026
# COMPLETED: 08 Sep 2026
# LAST MODIFIED: 08 Sep 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 0. Purpose ----
# ______________________________________________________________________________

# after BLASTing both the general grass and unidentified ASVs,
# we still have some that returned no matches. 
# Here, we'll separate those out from the ASVs we identified
# for another BLAST round

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(seqinr)

# ______________________________________________________________________________
# 2. Read data ----
# ______________________________________________________________________________

blast.grasses <- read.csv("BLAST/blast_grasses.csv")
blast.uni <- read.csv("BLAST/blast_unidentified.csv")

# ______________________________________________________________________________
# 3. Clean ----
# ______________________________________________________________________________

# split identified
blast.identified <- bind_rows(
  
  blast.grasses |> filter(group == "identified2"),
  blast.uni |> filter(group == "identified2")
  
)

# split ASVs for further BLASTing
blast.next <- bind_rows(
  
  blast.grasses |> filter(group == ""),
  blast.uni |> filter(group == "")
  
)

# ______________________________________________________________________________
# 4. Save to file ----
# ______________________________________________________________________________

write.csv(blast.identified, "data_cleaned/asvs_identified2.csv", row.names = F)

# ______________________________________________________________________________
# 5. FASTA file ----
# ______________________________________________________________________________

asvs <- read.csv("data_raw/asvs.csv") |>

  # unique ASV column
  mutate(ASV.unq = paste0(ASV, "_", batchID)) |>
  
  dplyr::select(sequence, ASV.unq) |>
  
  # right join
  right_join(blast.next)

# write
write.fasta(
  
  as.list(asvs$sequence),
  asvs$ASV.unq,
  "data_cleaned/FASTA/next.fasta",
  as.string = T
  
)

# ______________________________________________________________________________
# 6. After round 2 ----

# was able to identify some more ASVs by changing the BLAST parameters. 
# let's incorporate the new ASVs and separate out the "bad" ones from before
# to see if we can find more reasonable reads.

# these included Quercus and Betula; really shouldn't be in the hare zone

# ______________________________________________________________________________

# after next BLAST
blast.next <- read.csv("BLAST/blast_next.csv")

# split identified
blast.next.identified <- blast.next |> filter(group == "identified3")

# split unidentified
blast.next.uni <- blast.next |> filter(group == "")

# save to file
write.csv(blast.next.identified, "data_cleaned/asvs_identified3.csv", row.names = F)
write.csv(blast.next.uni, "data_cleaned/asvs_unidentified2.csv", row.names = F)

# "bad" reads
blast.uni <- read.csv("BLAST/blast_unidentified.csv")

blast.uni.bad <- blast.uni |> filter(group == "bad")

# FASTA
asvs <- read.csv("data_raw/asvs.csv") |>
  
  # unique ASV column
  mutate(ASV.unq = paste0(ASV, "_", batchID)) |>
  
  dplyr::select(sequence, ASV.unq) |>
  
  # right join
  right_join(blast.uni.bad)

# write
write.fasta(
  
  as.list(asvs$sequence),
  asvs$ASV.unq,
  "data_cleaned/FASTA/uni_bad.fasta",
  as.string = T
  
)

# ______________________________________________________________________________
# 6. After round 3 ----

# all "bad" IDs stayed the same, except for one (ASV_185_2, was Shepherdia)

# ______________________________________________________________________________

blast.bad <- read.csv("BLAST/blast_bad.csv")

blast.bad.id <- blast.bad |> filter(group == "identified4")

# save to file
write.csv(blast.bad.id, "data_cleaned/asvs_identified4.csv", row.names = F)
