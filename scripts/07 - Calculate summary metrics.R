# PROJECT: Diet
# SCRIPT: 07 - Calculate summary metrics
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 15 Sep 2026
# COMPLETED: 
# LAST MODIFIED: 15 Sep 2026
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
    
  )

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

# ______________________________________________________________________________
# 8. Wrapper function to calculate all ----
# ______________________________________________________________________________

calc_metabar_metrics <- function (.reads) {
  
  suppressMessages({
    
    foo <- calc_foo(.reads)
    poo <- calc_poo(.reads)
    wpoo <- calc_wpoo(.reads)
    rra <- calc_rra(.reads)
    
    # join
    all.metrics <- foo |> dplyr::select(final.taxon, pfoo) |>
      
      left_join(poo |> dplyr::select(final.taxon, poo)) |>
      
      left_join(wpoo |> dplyr::select(final.taxon, wpoo)) |>
      
      left_join(rra |> dplyr::select(final.taxon, rra))
    
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
  
  taxa.ranks <- data.frame(
    
    final.taxon = .metrics$final.taxon,
    rank.pfoo = rank(1 / .metrics$pfoo),
    rank.poo = rank(1 / .metrics$poo),
    rank.wpoo = rank(1 / .metrics$wpoo),
    rank.rra = rank(1 / .metrics$rra)
    
  ) |>
    
    # arrange by mean rank
    mutate(mean.rank = (rank.pfoo + rank.poo + rank.wpoo + rank.rra) / 4) |>
    
    arrange(mean.rank)
  
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

# NEXT:
  # functional groups
  # bootstrapping functions