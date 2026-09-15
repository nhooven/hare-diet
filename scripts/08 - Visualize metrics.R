# PROJECT: Diet
# SCRIPT: 08 - Visualize metrics
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 15 Sep 2026
# COMPLETED: 15 Sep 2026
# LAST MODIFIED: 15 Sep 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

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
  
  ungroup()

# join in
metrics.off <- metrics.off |> left_join(taxa.info)
metrics.on <- metrics.on |> left_join(taxa.info)

# ______________________________________________________________________________
# 3. Plot function ----
# ______________________________________________________________________________

plot_metric <- function (.metrics,
                         .which) {
  
  # keep correct column
  .metrics <- .metrics |> dplyr::select(final.taxon,
                                        all_of(.which),
                                        cat1,
                                        cat2) |>
    
    rename(metric = .which)
  
  # correct x-axis title
  x.title <- case_when(.which == "pfoo" ~ "% frequency of occurrence",
                       .which == "poo" ~ "Proportion of occurrence",
                       .which == "wpoo" ~ "Split-sample frequence of occurrence",
                       .which == "rra" ~ "Relative read abundance")
  
  # plot
  ggplot(data = .metrics) +
    
    theme_classic() +
    
    geom_col(aes(x = metric,
                 y = reorder(final.taxon, metric),
                 fill = cat1),
             linewidth = 0.2) +
    
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(color = "black",
                                     size = 8),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 10),
          legend.position = c(0.7, 0.4),
          legend.title = element_blank(),
          plot.margin = margin(t = 5, r = -1.5, b = 5, l = 5)) +
    
    scale_fill_viridis_d(option = "inferno") +
    
    scale_x_continuous(expand = c(0, 0)) +
    
    xlab(x.title)
  
}

# ______________________________________________________________________________
# 4. Plot ----
# ______________________________________________________________________________
# 4a. Snow-off ----

# 800 x 600

# ______________________________________________________________________________

plot_metric(metrics.off, "pfoo")
plot_metric(metrics.off, "poo")
plot_metric(metrics.off, "wpoo")
plot_metric(metrics.off, "rra")

# ______________________________________________________________________________
# 4b. Snow-off ----

# 800 x 400

# ______________________________________________________________________________

plot_metric(metrics.on, "pfoo")
plot_metric(metrics.on, "poo")
plot_metric(metrics.on, "wpoo")
plot_metric(metrics.on, "rra")
