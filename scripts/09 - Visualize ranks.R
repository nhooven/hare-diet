# PROJECT: Diet
# SCRIPT: 07 - Calculate summary metrics
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
    
    # factor levels
    mutate(name = factor(name,
                         levels = c("rank.pfoo",
                                    "rank.poo",
                                    "rank.wpoo",
                                    "rank.rra"),
                         labels = c("%FOO",
                                    "POO",
                                    "SSFOO",
                                    "RRA")))
  
  return(df.for.plot)
  
}

# use
ranks.off.for.plot <- prep_for_plot(ranks.off)
ranks.on.for.plot <- prep_for_plot(ranks.on)

# ______________________________________________________________________________
# 4. Plot ----
# ______________________________________________________________________________
# 4a. Function ----
# ______________________________________________________________________________

plot_ranks <- function(.df.for.plot,
                       .legend = T) {
  
  out.plot <- ggplot(data = .df.for.plot) +
    
    theme_classic() +
    
    geom_point(aes(x = value,
                   y = reorder(final.taxon, 1 / mean.rank),
                   shape = name,
                   size = name),
               fill = "white") +        # different color by functional group?
    
    theme(axis.text.y = element_text(size = 5),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_line(color = "gray95"),
          legend.position = c(0.2, 0.8),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "white",
                                           color = "black"),
          axis.ticks.x = element_blank()) +
    
    # legend
    guides(fill = "none",
           size = guide_legend(override.aes = list(size = 2.5))) +
    
    # axis
    scale_x_reverse(breaks = c(quantile(.df.for.plot$value, prob = 0.25),
                               quantile(.df.for.plot$value, prob = 0.75)),
                    labels = c("higher >",
                               "< lower")) +
    xlab("Importance rank") +
    coord_cartesian(xlim = c(2, max(.df.for.plot$value) - 2)) +
    
    # point shapes and sizes
    scale_shape_manual(values = c(0, 3, 4, 5)) +
    scale_size_manual(values = c(1.4, 1, 1.2, 1.2))
  
  if (.legend == F) {
    
    out.plot <- out.plot + theme(legend.position = "none")
    
  }
  
  return(out.plot)
  
}

# ______________________________________________________________________________
# 4b. Plot ----
# ______________________________________________________________________________

plot_ranks(ranks.off.for.plot)

# 800 x 600

plot_ranks(ranks.on.for.plot, .legend = F)

# 800 x 400