# PROJECT: Diet
# SCRIPT: 07 - Calculate summary metrics
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 15 Sep 2026
# COMPLETED: 15 Sep 2026
# LAST MODIFIED: 30 Sep 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(cowplot)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

ranks.off <- read.csv("data_cleaned/summaries/ranks_off.csv")
ranks.on <- read.csv("data_cleaned/summaries/ranks_on.csv")

# taxon information
taxa.info <- read.csv("data_cleaned/reads_taxa.csv") |>
  
  dplyr::select(final.taxon,
                cat1,
                cat2) |>
  
  group_by(final.taxon) |>
  
  slice(1) |>
  
  ungroup()

# scientific name set for labels
nonscientific.set <- taxa.info$final.taxon[c(15, 16, 21, 43, 59, 68, 76, 96, 97)]
scientific.set <- taxa.info$final.taxon[which(taxa.info$final.taxon %notin% nonscientific.set)]

sci_labels <- function (x) {
  
  lapply(x , function (lab) {
    
    if (lab %in% scientific.set) bquote(italic(.(lab))) else lab
    
  })
  
}

# ______________________________________________________________________________
# 3. Prepare for plotting ----
# ______________________________________________________________________________

# function
prep_for_plot <- function (.ranks) {
  
  df.for.plot <- .ranks |>
    
    # join in taxa.info
    left_join(taxa.info) |>
    
    # pivot
    pivot_longer(cols = c(rank.pfoo,
                          rank.poo,
                          rank.wpoo,
                          rank.rra)) |>
    
    # factor levels and functional groups
    mutate(name = factor(name,
                         levels = c("rank.pfoo",
                                    "rank.poo",
                                    "rank.wpoo",
                                    "rank.rra"),
                         labels = c("%FOO",
                                    "POO",
                                    "SSFOO",
                                    "RRA")),
           
           func = factor(cat1,
                         levels = c("conifer",
                                    "woody broadleaf",
                                    "sub-shrub",
                                    "forb", 
                                    "graminoid",
                                    "unknown")))
  
  return(df.for.plot)
  
}

# use
ranks.off.for.plot <- prep_for_plot(ranks.off)
ranks.on.for.plot <- prep_for_plot(ranks.on)

# colors
cat1.colors <- c("darkgreen", "darkorange3", "brown",
                 "green3", "lightgreen",
                 "gray")

# ______________________________________________________________________________
# 4. Plot ----
# ______________________________________________________________________________
# 4a. Function ----
# ______________________________________________________________________________

plot_ranks <- function(.df.for.plot,
                       .n.taxa = 30,
                       .legend = c("metric", "func")) {
  
  # only use the top n.taxa, which would be 4 metrics-worth
  .df.for.plot.1 <- .df.for.plot |> slice(1:(.n.taxa * 4))
  
  out.plot <- ggplot(data = .df.for.plot.1) +
    
    theme_classic() +
    
    geom_point(aes(x = value,
                   y = reorder(final.taxon, 1 / mean.rank),
                   shape = name,
                   size = name,
                   color = func),
               fill = "white") +        # different color by functional group?
    
    theme(axis.text.y = element_text(size = 5),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_line(color = "gray95"),
          legend.position = c(0.8, 0.2),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "white",
                                           color = "black"),
          legend.margin = margin(0.3, 4.5, 0.3, 0.5),
          legend.key.spacing.y = unit(-5, "pt"),
          legend.text = element_text(size = 8,
                                     margin = margin(-0.5, 0.5, -0.5, 0.5)),
          axis.ticks.x = element_blank()) +
    
    # axis
    scale_x_reverse(breaks = c(max(.df.for.plot.1$value - 3),
                               min(.df.for.plot.1$value + 3)),
                    labels = c("low", "high")) +
    xlab("Importance rank") +
    coord_cartesian(xlim = c(2, max(.df.for.plot.1$value) - 2)) +
    
    scale_y_discrete(labels = sci_labels) +
    
    # point shapes and sizes
    scale_shape_manual(values = c(0, 3, 4, 5)) +
    scale_size_manual(values = c(1.4, 1, 1.2, 1.2)) +
    scale_color_manual(values = cat1.colors)
  
  if (.legend == "metric") {
    
    out.plot <- out.plot + 
      
      # legend
      guides(fill = "none",
             size = guide_legend(override.aes = list(size = 2))) +
      
      scale_color_manual(values = cat1.colors,
                         guide = "none")
    
  }
  
  if (.legend == "func") {
    
    out.plot <- out.plot + 
    
    # legend
    guides(color = guide_legend(override.aes = list(size = 2))) +
      
      scale_size_manual(values = c(1.4, 1, 1.2, 1.2),
                        guide = "none") +
      scale_shape_manual(values = c(0, 3, 4, 5),
                         guide = "none")
    
  }
  
  return(out.plot)
  
}

# ______________________________________________________________________________
# 4b. Plot ----

# 850 x 380

# ______________________________________________________________________________

plot_grid(plot_ranks(ranks.off.for.plot, 30, "metric"),
          plot_ranks(ranks.on.for.plot, 30, "func"),
          labels = "auto",
          nrow = 1)
