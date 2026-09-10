# PROJECT: Diet
# SCRIPT: 06 - Data exploration
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 09 Sep 2026
# COMPLETED: 
# LAST MODIFIED: 10 Sep 2026
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
# 3. Attribute sample data ----
# ______________________________________________________________________________

reads.1 <- reads.taxa |>
  
  # add in Ear.tag, Treatment, Season, Sex, Site
  left_join(
    
    samples.lookup |> dplyr::select(Final.sample.ID, 
                                    Ear.tag, Treatment, Season, Sex, Site)
    
  )

# ______________________________________________________________________________
# 4. Split by season ----
# ______________________________________________________________________________

reads.off <- reads.1 |> filter(Season == "Off")
reads.on <- reads.1 |> filter(Season == "On")

# ______________________________________________________________________________
# 5. Remove taxa <=1% of reads in a given sample ----

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
reads.off.1 <- remove_rare_taxa(reads.off)
reads.on.1 <- remove_rare_taxa(reads.on)

# ______________________________________________________________________________
# 5. Functions to calculate metrics ----
# ______________________________________________________________________________
# 5a. Frequency of occurrence (FOO) ----

# Deagle et al 2019:
  # the number of samples that contain a given food item,
  # most often expressed as a per cent

# ______________________________________________________________________________

calc_foo <- function (.reads) {
  
  # unique food items
  unique.items <- unique(.reads$final.taxon)
  n.items <- length(unique.items)
  
  # number of samples
  unique.samples <- unique(.reads$Final.sample.ID)
  n.samples <- length(unique.samples)
  
  # helper function - check how many samples each taxon is in
  check_foo <- function (.taxon) {
    
    length(unique(.reads$Final.sample.ID[.reads$final.taxon == .taxon]))
    
  }
  
  # list and lapply
  unique.items.list <- as.list(unique.items)
  
  taxa.foo <- data.frame(
    
    final.taxon = unique.items,
    foo = do.call(rbind, lapply(unique.items.list, check_foo))
    
  ) |>
    
    mutate(pfoo = foo / n.samples) |>
    
    # arrange by pfoo
    arrange(desc(pfoo))
  
  return(taxa.foo)
  
}

# ______________________________________________________________________________
# 5b. Percent of occurrence (POO) ----

# Deagle et al 2019:
  # simply pFOO rescaled so that the sum across all food items is 100%

# ______________________________________________________________________________

calc_poo <- function (.reads) {
  
  # call calc_foo
  taxa.foo <- calc_foo(.reads)
  
  # rescale
  taxa.poo <- taxa.foo |>
    
    dplyr::select(final.taxon, foo) |>
    
    mutate(poo = foo / sum(foo))
  
  return(taxa.poo)
  
}

# ______________________________________________________________________________
# 5c. Weighted percent of occurrence (wPOO) ----

# Deagle et al 2019:
  # similar to POO, but rather than giving equal weight to all occurrences, 
  # this metric weights each occurrence according to the number of food items in the sample 
  # (e.g., if a sample contains five food items, each will be given weight 1/5)

# ______________________________________________________________________________

calc_wpoo <- function (.reads) {
  
  # unique food items
  unique.items <- unique(.reads$final.taxon)
  n.items <- length(unique.items)
  
  # number of samples
  unique.samples <- unique(.reads$Final.sample.ID)
  n.samples <- length(unique.samples)
  
  # helper function - return the total items by sample
  return_items_in_sample <- function (.sample) {
    
    length(unique(.reads$final.taxon[.reads$Final.sample.ID == .sample]))
    
  }
  
  # lapply by sample
  unique.samples.list <- as.list(unique.samples)
  
  sample.n.items <- data.frame(
    
    Final.sample.ID = unique.samples,
    sample.n.items = do.call(rbind, lapply(unique.samples.list, return_items_in_sample))
    
  )
  
  # expand with taxa
  taxa.wpoo <- sample.n.items |>
    
    left_join(
      
      .reads |> dplyr::select(Final.sample.ID,
                              final.taxon)
      
    ) |>
    
    # create weights
    mutate(weight = 1 / sample.n.items) |>
    
    # sum weights by taxon
    group_by(final.taxon) |>
    summarize(wsum = sum(weight))|>
    ungroup() |>
    
    # wPOO - divide by total samples
    mutate(wpoo = wsum / n.samples) |>
    
    # arrange
    arrange(desc(wpoo))
  
  return(taxa.wpoo)
  
}

# ______________________________________________________________________________
# 5d. Relative read abundance (RRA) ----

# Deagle et al 2019:
  # item reads / total reads in sample, summed over all samples, divided by total samples

# ______________________________________________________________________________

calc_rra <- function (.reads) {
  
  # unique food items
  unique.items <- unique(.reads$final.taxon)
  n.items <- length(unique.items)
  
  # number of samples
  unique.samples <- unique(.reads$Final.sample.ID)
  n.samples <- length(unique.samples)
  
  # total reads per sample
  total.reads.sample <- .reads |>
    
    group_by(Final.sample.ID) |>
    summarize(total.reads = sum(reads)) |>
    ungroup()
  
  # RRA
  taxa.rra <- .reads |> 
      
    # total reads per taxon, per sample
    group_by(final.taxon, Final.sample.ID) |>
    summarize(item.reads = sum(reads)) |>
    ungroup() |>
    
    # join in
    left_join(total.reads.sample) |>
    
    # divide
    mutate(rel.reads = item.reads / total.reads) |>
    
    # sum by taxon
    group_by(final.taxon) |>
    summarize(sum.rel.reads = sum(rel.reads)) |>
    ungroup() |>
    
    # divide by total samples
    mutate(rra = sum.rel.reads / n.samples) |>
    
    # arrange
    arrange(desc(rra))
  
  return(taxa.rra)
  
}

# ______________________________________________________________________________
# 5e. Wrapper function - All metrics ----
# ______________________________________________________________________________

calc_metabar_metrics <- function (.reads) {
  
  suppressMessages({
    
    foo <- calc_foo(.reads)
    poo <- calc_poo(.reads)
    wpoo <- calc_wpoo(.reads)
    rra <- calc_rra(.reads)
    
    # join
    all.metrics <- foo |>
      
      left_join(poo |> dplyr::select(-foo)) |>
      
      left_join(wpoo |> dplyr::select(-wsum)) |>
      
      left_join(rra |> dplyr::select(-sum.rel.reads))
    
  })
  
  return(all.metrics)
  
}

# ______________________________________________________________________________
# 6. Calculate metrics ----
# ______________________________________________________________________________
# 6a. By season ----
# ______________________________________________________________________________

metrics.off <- calc_metabar_metrics(reads.off.1)
metrics.on <- calc_metabar_metrics(reads.on.1)

# ______________________________________________________________________________
# 7. Bar charts ----

library(cowplot)

# ______________________________________________________________________________
# 7a. By season ----
# ______________________________________________________________________________

# snow-off
bar.off.pfoo <- ggplot(data = metrics.off) +
  
  theme_classic() +
  
  geom_col(aes(x = pfoo,
               y = reorder(final.taxon, pfoo),
               fill = pfoo),
           color = "darkgray",
           linewidth = 0.2) +
  
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        legend.position = "none",
        plot.margin = margin(t = 5, r = -1.5, b = 5, l = 5)) +
  
  xlab("% frequency of occurrence") +
  
  scale_fill_viridis_c(option = "inferno") +
  
  # reverse the x axis
  scale_x_reverse(breaks = c(0.2, 0.4, 0.6, 0.8),
                  labels = c(20, 40, 60, 80)) +
  
  # move axis labels to other side
  scale_y_discrete(position = "right") +
  
  # axis limits
  coord_cartesian(xlim = c(0.047, max(metrics.off$pfoo)))

# wPOO
bar.off.wpoo <- ggplot(data = metrics.off) +
  
  theme_classic() +
  
  geom_col(aes(x = wpoo,
               y = reorder(final.taxon, pfoo),
               fill = wpoo),
           color = "darkgray",
           linewidth = 0.2) +
  
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(color = "black",  # must use text in this one
                                   size = 5,
                                   hjust = 0.5),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        legend.position = "none",
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        plot.background = element_rect(fill = NA,
                                       color = NA)) +
  
  xlab("Weighted percent of occurrence") +
  
  scale_fill_viridis_c(option = "inferno") +
  
  # axis limits
  coord_cartesian(xlim = c(0.006, max(metrics.off$wpoo))) +
  scale_x_continuous(breaks = c(0.02, 0.04, 0.06, 0.08, 0.1, 0.12),
                     labels = c(2, 4, 6, 8, 10, 12))

# plot together
plot_grid(bar.off.pfoo, bar.off.wpoo,
          nrow = 1,
          rel_widths = c(1, 1.4))

# 650 x 550

# snow-on
bar.on.pfoo <- ggplot(data = metrics.on) +
  
  theme_classic() +
  
  geom_col(aes(x = pfoo,
               y = reorder(final.taxon, pfoo),
               fill = pfoo),
           color = "darkgray",
           linewidth = 0.2) +
  
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        legend.position = "none",
        plot.margin = margin(t = 5, r = -1.5, b = 5, l = 5)) +
  
  xlab("% frequency of occurrence") +
  
  scale_fill_viridis_c(option = "mako") +
  
  # reverse the x axis
  scale_x_reverse(breaks = c(0.2, 0.4, 0.6, 0.8),
                  labels = c(20, 40, 60, 80)) +
  
  # move axis labels to other side
  scale_y_discrete(position = "right") +
  
  # axis limits
  coord_cartesian(xlim = c(0.047, max(metrics.on$pfoo)))

# wPOO
bar.on.wpoo <- ggplot(data = metrics.on) +
  
  theme_classic() +
  
  geom_col(aes(x = wpoo,
               y = reorder(final.taxon, pfoo),
               fill = wpoo),
           color = "darkgray",
           linewidth = 0.2) +
  
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(color = "black",  # must use text in this one
                                   size = 5,
                                   hjust = 0.5),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        legend.position = "none",
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        plot.background = element_rect(fill = NA,
                                       color = NA)) +
  
  xlab("Weighted percent of occurrence") +
  
  scale_fill_viridis_c(option = "mako") +
  
  # axis limits
  coord_cartesian(xlim = c(0.046, max(metrics.on$wpoo))) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0),
                     labels = c(20, 40, 60, 80, 100))

# plot together
plot_grid(bar.on.pfoo, bar.on.wpoo,
          nrow = 1,
          rel_widths = c(1, 1.4))

# 650 x 367
