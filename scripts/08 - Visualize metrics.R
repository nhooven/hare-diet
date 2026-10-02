# PROJECT: Diet
# SCRIPT: 08 - Visualize metrics
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 15 Sep 2026
# COMPLETED: 15 Sep 2026
# LAST MODIFIED: 02 Oct 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(cowplot)
library(mefa4)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

metrics.off <- read.csv("data_cleaned/summaries/metrics_off.csv")
metrics.on <- read.csv("data_cleaned/summaries/metrics_on.csv")

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

# join in
metrics.off <- metrics.off |> left_join(taxa.info)
metrics.on <- metrics.on |> left_join(taxa.info)

# scientific name set for labels
nonscientific.set <- taxa.info$final.taxon[c(15, 16, 21, 43, 59, 68, 76, 96, 97)]
scientific.set <- taxa.info$final.taxon[which(taxa.info$final.taxon %notin% nonscientific.set)]

sci_labels <- function (x) {
  
  lapply(x , function (lab) {
    
    if (lab %in% scientific.set) bquote(italic(.(lab))) else lab
    
  })
  
}

# ______________________________________________________________________________
# 3. Plot function ----

# colors
cat1.colors <- c("darkgreen", "darkorange3", "brown",
                 "green3", "lightgreen",
                 "gray")

# ______________________________________________________________________________

plot_metric <- function (.metrics,
                         .which,
                         .n.taxa = 35,
                         .legend = T) {
  
  # keep correct column
  .metrics <- .metrics |> filter(metric == .which) |>
    
    dplyr::select(final.taxon,
                  value,
                  lci,
                  uci,
                  cat1,
                  cat2,
                  func) |>
    
    mutate(metric = .which) |>
    
    # keep only the first n
    slice(1:.n.taxa)
  
  # correct x-axis title
  x.title <- case_when(.which == "pfoo" ~ "% frequency of occurrence",
                       .which == "poo" ~ "Proportion of occurrence",
                       .which == "wpoo" ~ "Split-sample frequency of occurrence",
                       .which == "rra" ~ "Relative read abundance")
  
  # plot
  out.plot <- ggplot(data = .metrics) +
    
    theme_classic() +
    
    geom_col(aes(x = value,
                 y = reorder(final.taxon, value),
                 fill = func),
             linewidth = 0.2) +
    
    geom_errorbar(aes(x = value,
                      y = reorder(final.taxon, value),
                      xmin = lci,
                      xmax = uci),
                  height = 0,
                  color = "darkgray") +
    
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(color = "black",
                                     size = 8),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 9),
          legend.position = c(0.7, 0.35),
          legend.title = element_blank(),
          legend.text = element_text(size = 8,
                                     vjust = 1),
          legend.key.size = unit(0.35, "cm"),
          plot.margin = margin(t = 1, r = 10, b = 1, l = 7)) +
    
    scale_fill_manual(values = cat1.colors) +
    
    scale_y_discrete(labels = sci_labels) +
    
    scale_x_continuous(expand = c(0, 0)) +
    
    xlab(x.title) 
  
  if (.legend == F) {
    
    out.plot <- out.plot + theme(legend.position = "none")
    
  }
  
  return(out.plot)
  
}

# ______________________________________________________________________________
# 4. Plot ----

# 620 x 480

# ______________________________________________________________________________
# 4a. Snow-off ----
# ______________________________________________________________________________

plot_grid(plot_metric(metrics.off, "pfoo", 30, T),
          plot_metric(metrics.off, "poo", 30, F),
          plot_metric(metrics.off, "wpoo", 30, F),
          plot_metric(metrics.off, "rra", 30, F),
          nrow = 2)

# ______________________________________________________________________________
# 4b. Snow-on ----
# ______________________________________________________________________________

plot_grid(plot_metric(metrics.on, "pfoo", 30, T),
          plot_metric(metrics.on, "poo", 30, F),
          plot_metric(metrics.on, "wpoo", 30, F),
          plot_metric(metrics.on, "rra", 30, F),
          nrow = 2)
