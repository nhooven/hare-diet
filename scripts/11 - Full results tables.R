# PROJECT: Diet
# SCRIPT: 11 - Full results tables
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 02 Oct 2026
# COMPLETED: 02 Oct 2026
# LAST MODIFIED: 02 Oct 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 0. Purpose ----
# ______________________________________________________________________________

# this script will create complete results tables (rather than figures with only top
# 30-or-so taxa), along with the 
  # total reads for that taxon
  # total samples that include it
  # calculated metric, SE, and rank

# sorted by mean rank?

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# metrics
metrics.off <- read.csv("data_cleaned/summaries/metrics_off.csv")
metrics.on <- read.csv("data_cleaned/summaries/metrics_on.csv")

# ranks
ranks.off <- read.csv("data_cleaned/summaries/ranks_off.csv")
ranks.on <- read.csv("data_cleaned/summaries/ranks_on.csv")

# cleaned reads
reads.clean <- read.csv("data_cleaned/reads_cleaned.csv")

# foo function
source("functions/calc_foo.R")

# taxon information
taxa.info <- read.csv("data_cleaned/reads_taxa.csv") |>
  
  dplyr::select(final.taxon,
                cat1,
                cat2) |>
  
  group_by(final.taxon) |>
  
  slice(1) |>
  
  ungroup() |>
  
  # functional groups
  mutate(func = factor(cat1,
                       levels = c("conifer",
                                  "woody broadleaf",
                                  "sub-shrub",
                                  "forb", 
                                  "graminoid",
                                  "unknown")))

# ______________________________________________________________________________
# 3. Function - create table ----
# ______________________________________________________________________________

create_table <- function (.season = c("off", "on")) {
  
  # correct datasets
  if (.season == "off") { 
    
    .metrics <- metrics.off
    .ranks <- ranks.off 
    .reads <- reads.clean |> filter(Season == "Off")
    
  } else {
      
    .metrics <- metrics.on
    .ranks <- ranks.on
    .reads <- reads.clean |> filter(Season == "On")
    
  }
  
  # total reads by taxon
  total.taxon.reads <- .reads |>
    
    group_by(final.taxon) |>
    
    summarize(total.reads = sum(reads))
  
  # sample frequency of occurrence
  total.taxon.samples <- calc_foo(.reads) |>
    
    dplyr::select(final.taxon, foo)
  
  # create table
  full.table <- .metrics |>
    
    # remove CIs
    dplyr::select(-c(lci, uci)) |>
    
    # pivot_wider
    pivot_wider(names_from = metric,
                values_from = c(value, se)) |>
    
    # round to 3 (?) digits
    mutate(across(value_pfoo:se_rra, \(x) round(x, digits = 3))) |>
    
    # join
      # ranks
      left_join(.ranks) |>
    
      # functional group
      left_join(taxa.info |> select(final.taxon, cat1)) |>
    
      # reads
      left_join(total.taxon.reads) |>
      
      # n samples (foo)
      left_join(total.taxon.samples) |>
    
    # arrange by mean rank
    arrange(mean.rank) |>
    
    # and put columns in correct order
    dplyr::select(final.taxon,
                  cat1,
                  foo,
                  total.reads,
                  value_pfoo, se_pfoo, rank.pfoo,
                  value_poo, se_poo, rank.poo,
                  value_wpoo, se_wpoo, rank.wpoo,
                  value_rra, se_rra, rank.rra)
  
  return(full.table)
  
}

# ______________________________________________________________________________
# 4. Create and write tables ----
# ______________________________________________________________________________

write.table(create_table("off"), "clipboard", sep = "\t")
write.table(create_table("on"), "clipboard", sep = "\t")            
