# PROJECT: Diet
# SCRIPT: 07 - Calculate summary metrics
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 15 Sep 2026
# COMPLETED: 02 Oct 2026
# LAST MODIFIED: 02 Oct 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

samples.lookup <- read.csv("data_cleaned/samples_lookup.csv")
reads.taxa <- read.csv("data_cleaned/reads_taxa.csv")

# ______________________________________________________________________________
# 3. Keep samples within cutoff ----
# ______________________________________________________________________________

samples.lookup.cut1 <- samples.lookup |> 
  filter(include.reasonable == "Y") |> 
  dplyr::select(-c(include.reasonable, include.stringent))

samples.lookup.cut2 <- samples.lookup |> 
  filter(include.stringent == "Y") |> 
  dplyr::select(-c(include.reasonable, include.stringent))

# ______________________________________________________________________________
# 4. Attribute sample data ----

# we'll start with the "reasonable" cutoff

# ______________________________________________________________________________

reads.1 <- reads.taxa |>
  
  # add in Ear.tag, Treatment, Season, Sex, Site
  left_join(
    
    samples.lookup.cut1 |> dplyr::select(Final.sample.ID, 
                                         Ear.tag, Treatment, Season, Sex, Site)
    
  ) |>
  
  # keep only complete cases
  drop_na(Ear.tag)

# ______________________________________________________________________________
# 5. Remove taxa <= 1% of reads in a given sample ----

# Deagle et al 2019

# ______________________________________________________________________________

# function
remove_rare_taxa <- function (.reads) {
  
  # number of samples
  unique.samples <- unique(.reads$Final.sample.ID)
  n.samples <- length(unique.samples)
  
  # helper function - calculate percent of reads by taxa, within sample
  calc_perc_reads <- function (.sample) {
    
    sample.reads <- .reads |> filter(Final.sample.ID == .sample)
    
    sample.reads.1 <- sample.reads |>
      
      group_by(final.taxon) |>
      summarize(reads.by.taxon = sum(reads)) |>
      ungroup() |>
      
      mutate(perc.total.reads = reads.by.taxon / sum(reads.by.taxon)) |>
      
      # filter only those with >1% of reads
      filter(perc.total.reads > 0.01)
    
    # return sample reads with rare taxa removed
    sample.reads.out <- sample.reads |>
      
      filter(final.taxon %in% sample.reads.1$final.taxon)
    
  }
  
  # lapply
  unique.samples.list <- as.list(unique.samples)
  
  return(do.call(rbind, lapply(unique.samples.list, calc_perc_reads)))
  
}

# use
reads.2 <- remove_rare_taxa(reads.1)

# ______________________________________________________________________________
# 6. Split dataset ----
# ______________________________________________________________________________

# season
reads.off <- reads.2 |> filter(Season == "Off")
reads.on <- reads.2 |> filter(Season == "On")

# ______________________________________________________________________________
# 7. Read metric functions ----
# ______________________________________________________________________________

source("functions/calc_foo.R")
source("functions/calc_poo.R")
source("functions/calc_wpoo.R")
source("functions/calc_rra.R")

source("functions/boot_metric.R")

# ______________________________________________________________________________
# 8. Wrapper function to calculate all ----
# ______________________________________________________________________________

calc_metabar_metrics <- function (.reads) {
  
  suppressMessages({
    
    foo <- boot_metric(.reads, .iter = 1000, "foo")
    poo <- boot_metric(.reads, .iter = 1000, "poo")
    wpoo <- boot_metric(.reads, .iter = 1000, "wpoo")
    rra <- boot_metric(.reads, .iter = 1000, "rra")
    
    # join (create metric column too)
    all.metrics <- foo |> dplyr::select(final.taxon, pfoo, se, lci, uci) |> 
      
      rename(value = pfoo) |>
      mutate(metric = "pfoo") |>
      
      bind_rows(poo |> dplyr::select(final.taxon, poo, se, lci, uci) |> 
                  
                  rename(value = poo) |>
                  mutate(metric = "poo")) |>
      
      bind_rows(wpoo |> dplyr::select(final.taxon, wpoo, se, lci, uci) |> 
                  
                  rename(value = wpoo) |>
                  mutate(metric = "wpoo")) |>
      
      bind_rows(rra |> dplyr::select(final.taxon, rra, se, lci, uci) |> 
                  
                  rename(value = rra) |>
                  mutate(metric = "rra"))
    
  })
  
  return(all.metrics)
  
}

# calculate metrics
metrics.off <- calc_metabar_metrics(reads.off)
metrics.on <- calc_metabar_metrics(reads.on)

# ______________________________________________________________________________
# 9. Calculate importance ranks ----
# ______________________________________________________________________________

# function
calc_ranks <- function(.metrics) {
  
  taxa.ranks <- .metrics |>
    
    # ranks only need the value
    dplyr::select(final.taxon, value, metric) |>
    
    # pivot
    pivot_wider(names_from = "metric",
                values_from = "value") |>
    
    # add ranks
    mutate(
      
      rank.pfoo = rank(1 / pfoo),
      rank.poo = rank(1 / poo),
      rank.wpoo = rank(1 / wpoo),
      rank.rra = rank(1 / rra)
      
    ) |>
    
    # arrange by mean rank
    mutate(mean.rank = (rank.pfoo + rank.poo + rank.wpoo + rank.rra) / 4) |>
    
    arrange(mean.rank) |>
    
    # keep only taxa and ranks
    dplyr::select(final.taxon, rank.pfoo:rank.rra, mean.rank)
  
  return(taxa.ranks)
  
}

# calculate
ranks.off <- calc_ranks(metrics.off)
ranks.on <- calc_ranks(metrics.on)

# ______________________________________________________________________________
# 10. Write to files ----
# ______________________________________________________________________________

# cleaned reads
write.csv(reads.2, "data_cleaned/reads_cleaned.csv", row.names = F)

# metrics
write.csv(metrics.off, "data_cleaned/summaries/metrics_off.csv", row.names = F)
write.csv(metrics.on, "data_cleaned/summaries/metrics_on.csv", row.names = F)

# ranks
write.csv(ranks.off, "data_cleaned/summaries/ranks_off.csv", row.names = F)
write.csv(ranks.on, "data_cleaned/summaries/ranks_on.csv", row.names = F)
