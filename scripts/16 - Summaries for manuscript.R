# PROJECT: Diet
# SCRIPT: 16 - Summaries for manuscript
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 08 Oct 2026
# COMPLETED: 08 Oct 2026
# LAST MODIFIED: 08 Oct 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Samples ----
# ______________________________________________________________________________

# how many samples were successfully sequenced?
# how many did we remove with the reasonable cutoff? 

# sample lookup table
samples.lookup <- read.csv("data_cleaned/samples_lookup.csv")

nrow(samples.lookup)  # total samples = 183
sum(samples.lookup$include.reasonable == "Y")   # retained samples = 178

# how many did we get reads from?
reads <- read.csv("data_cleaned/reads_cleaned.csv")

length(unique(reads$Final.sample.ID))   # apparently 178

# per season
nrow(samples.lookup |> filter(include.reasonable == "Y" & Season == "Off"))
nrow(samples.lookup |> filter(include.reasonable == "Y" & Season == "On"))

# ______________________________________________________________________________
# 3. Reads ----
# ______________________________________________________________________________

# total reads before and after filtering
# reads per sample

# all reads
reads.all <- read.csv("data_cleaned/reads_all.csv")

# total reads before filtering
sum(reads.all$reads)   # total reads = 3,911,158

# total reads after filtering
sum(reads$reads)       # total reads = 3,005,528

# reads per sample (with contaminants after threshold)
reads.per.sample.all <- reads.all |>
  
  group_by(Final.sample.ID) |>
  summarize(reads = sum(reads)) |>
  ungroup() |>
  
  filter(reads > 792)

mean(reads.per.sample.all$reads)     # mean = 21966.46
sd(reads.per.sample.all$reads)       # sd = 29104.9
range(reads.per.sample.all$reads)    # range = 1028 - 187386

# reads per sample (non-contaminant)
reads.per.sample <- reads |>
  
  group_by(Final.sample.ID) |>
  summarize(reads = sum(reads)) |>
  ungroup()
  
mean(reads.per.sample$reads)     # mean = 16884.99
sd(reads.per.sample$reads)       # sd = 23901.93
range(reads.per.sample$reads)    # range = 292 - 178248

# ______________________________________________________________________________
# 4. ASVs ----
# ______________________________________________________________________________

# total ASVs
# discarded and retained

asvs <- read.csv("data_raw/asvs.csv")
taxa <- read.csv("data_cleaned/all_taxa.csv")

length(unique(paste0(asvs$ASV, asvs$batchID)))  # total ASVs = 385
length(unique(taxa$ASV.unq))                    # retained ASVs = 333

# total unknown ASVs
nrow(taxa |> filter(final.taxon == "unknown"))  # total = 40

# meaning that we assigned 333 - 40 to at least a family
(333 - 40) / 333

# ______________________________________________________________________________
# 5. Diet items ----
# ______________________________________________________________________________

# total diet items
# diet items per sample
# how many to family, genus, species

length(unique(taxa$final.taxon))    # total taxa = 97 (96 + unknown)

# by season
reads.off <- reads |> filter(Season == "Off")
reads.on <- reads |> filter(Season == "On")

length(unique(reads.off$final.taxon))    # 60 during summer
length(unique(reads.on$final.taxon))     # 40 during winter

# by sample
items.sample.off <- reads.off |>
  
  group_by(Final.sample.ID) |>
  summarize(total.items = length(unique(final.taxon))) |>
  ungroup()

items.sample.on <- reads.on |>
  
  group_by(Final.sample.ID) |>
  summarize(total.items = length(unique(final.taxon))) |>
  ungroup()

# off
mean(items.sample.off$total.items)    # mean = 8.63
sd(items.sample.off$total.items)      # sd = 2.46
range(items.sample.off$total.items)   # range = 4 - 15

# on
mean(items.sample.on$total.items)    # mean = 5.88
sd(items.sample.on$total.items)      # sd = 2.27
range(items.sample.on$total.items)   # range = 2 - 18

# total to each rank
one.obs.per.taxon <- taxa |> group_by(final.taxon) |> slice(1) |> filter(final.taxon != "unknown")

# family 100%

# genus (-3 because of ambiguity)
nrow(one.obs.per.taxon |> filter(genus != "" & is.na(genus) == F))  # 87 - 93.75%

# species
nrow(one.obs.per.taxon |> filter(species != "" & is.na(species) == F))  # 56 - 58.3%
