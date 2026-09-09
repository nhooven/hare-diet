# PROJECT: Diet
# SCRIPT: 04 - Attribute reads to samples
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 09 Sep 2026
# COMPLETED: 09 Sep 2026
# LAST MODIFIED: 09 Sep 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 0. Metadata ----
# ______________________________________________________________________________

# samples_fieldID
  # my sample info sheet, including seasons, dates, sex, etc.

# samples_finalID
  # lookup table relating my fieldID with the Sacks Lab finalID

# reads_batchX
  # ASV reads per sample (batch 2 is transposed)

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# sample info
samples.fieldID <- read.csv("data_raw/samples_fieldID.csv")
samples.finalID <- read.csv("data_raw/samples_finalID.csv")

# reads
reads.batch1 <- read.csv("data_raw/reads_batch1.csv")
reads.batch2 <- read.csv("data_raw/reads_batch2.csv")
reads.batch3 <- read.csv("data_raw/reads_batch3.csv")

# ______________________________________________________________________________
# 3. Clean sample lookup ----
# ______________________________________________________________________________

# keep only the "new" replicates for batch3
# looks like the lab ran these twice - we'll keep the second one
which(samples.finalID$Final.sample.ID %in% c("S25-5592_R1",
                                             "S25-5593_R1",
                                             "S25-5594_R1",
                                             "S25-5595_R1",
                                             "S25-5596_R1",
                                             "S25-5597_R1",
                                             "S25-5598_R1",
                                             "S25-5599_R1")) -> which.dup

samples.finalID <- samples.finalID[-which.dup, ]

# keep consistent finalID format
samples.finalID <- samples.finalID |>
  
  mutate(Final.sample.ID = paste0(substr(Final.sample.ID,
                                         1, 3),
                                  "-",
                                  substr(Final.sample.ID,
                                         5, 8)))

# full lookup
str(samples.fieldID)

samples.lookup <- samples.fieldID |>
  
  # keep relevant columns
  dplyr::select(Batch,
                Sample,
                Period,
                Treatment,
                Season,
                Site,
                Ear.tag,
                Sex,
                Collection.date) |>
  
  # rename
  rename(FieldID = Sample) |>
  
  # coerce to date (why not)
  mutate(Collection.date = mdy(Collection.date)) |>
  
  # join in finalID
  left_join(samples.finalID)
  
# ______________________________________________________________________________
# 4. Clean reads ----
# ______________________________________________________________________________

# batch1
reads.batch1.1 <- reads.batch1 |>
  
  # unique ASV
  mutate(ASV.unq = paste0(ASV, "_", Batch)) |>
  
  # drop Batch and ASV (would complicate joining/pivoting)
  dplyr::select(-c(Batch, ASV)) |>
  
  # pivot longer
  pivot_longer(cols = 1:6,
               names_to = "Final.sample.ID",
               values_to = "reads") |>
  
  # correct name
  mutate(Final.sample.ID = paste0(substr(Final.sample.ID,
                                         1, 3),
                                  "-",
                                  substr(Final.sample.ID,
                                         5, 8))) |>
  
  # remove zero reads
  filter(reads > 0) |>
  
  # rearrange
  dplyr::select(Final.sample.ID,
                ASV.unq,
                reads)

# batch2
reads.batch2 <- reads.batch2 |> drop_na()

reads.batch2.mat <- reads.batch2 |>
  
  # remove batch and sampleID
  dplyr::select(-c(batch, sample_finalID)) |>
  
  # transpose
  t() |>
  
  # df
  as.data.frame()

# new row and column names
reads.batch2.1 <- cbind(data.frame(batch = 2),
                        colnames(reads.batch2[-c(1:2)]),
                        reads.batch2.mat)

colnames(reads.batch2.1) <- c("Batch", "ASV", reads.batch2$sample_finalID)

# clean
reads.batch2.1 <- reads.batch2.1 |>
  
  # unique ASV
  mutate(ASV.unq = paste0(ASV, "_", Batch)) |>
  
  # drop Batch and ASV (would complicate joining/pivoting)
  dplyr::select(-c(Batch, ASV)) |>
  
  # pivot longer
  pivot_longer(cols = 1:162,
               names_to = "Final.sample.ID",
               values_to = "reads") |>
  
  # correct name
  mutate(Final.sample.ID = paste0(substr(Final.sample.ID,
                                         1, 3),
                                  "-",
                                  substr(Final.sample.ID,
                                         5, 8))) |>
  
  # remove zero reads
  filter(reads > 0) |>
  
  # rearrange
  dplyr::select(Final.sample.ID,
                ASV.unq,
                reads)

# batch3
reads.batch3.1 <- reads.batch3 |>
  
  # unique ASV
  mutate(ASV.unq = paste0(ASV, "_", batch)) |>
  
  # drop Batch and ASV (would complicate joining/pivoting)
  dplyr::select(-c(batch, ASV)) |>
  
  # pivot longer
  pivot_longer(cols = 1:23,
               names_to = "Final.sample.ID",
               values_to = "reads") |>
  
  # correct name
  mutate(Final.sample.ID = paste0(substr(Final.sample.ID,
                                         1, 3),
                                  "-",
                                  substr(Final.sample.ID,
                                         5, 8))) |>
  
  # remove zero reads
  filter(reads > 0) |>
  
  # rearrange
  dplyr::select(Final.sample.ID,
                ASV.unq,
                reads)

# bind together
reads.all <- bind_rows(
  
  reads.batch1.1 |> mutate(batch = 1),
  reads.batch2.1 |> mutate(batch = 2),
  reads.batch3.1 |> mutate(batch = 3)
  
)

# ______________________________________________________________________________
# 5. Save to file ----
# ______________________________________________________________________________

write.csv(samples.lookup, "data_cleaned/samples_lookup.csv", row.names = F)
write.csv(reads.all, "data_cleaned/reads_all.csv", row.names = F)
