# PROJECT: Diet
# SCRIPT: 08 - Visualize metrics
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 15 Sep 2026
# COMPLETED: 15 Sep 2026
# LAST MODIFIED: 26 Oct 2026
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
nonscientific.set <- taxa.info$final.taxon[c(15, 16, 21, 42, 59, 68, 76, 96)]
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
  x.title <- case_when(.which == "pfoo" ~ "Frequency of occurrence",
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
          legend.position = c(0.7, 0.3),
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

# ______________________________________________________________________________
# 5. Pairwise correlation between metrics ----
# ______________________________________________________________________________

# clean
metric.corr.off <- metrics.off |>
  
  dplyr::select(final.taxon,
                value,
                metric) |>
  
  pivot_wider(values_from = value,
              names_from = metric)

metric.corr.on <- metrics.on |>
  
  dplyr::select(final.taxon,
                value,
                metric) |>
  
  pivot_wider(values_from = value,
              names_from = metric)

# Spearman's correlation
write.table(cor(metric.corr.off[ , 2:5], method = "spearman"), 
            "clipboard",
            sep = "\t")

write.table(cor(metric.corr.on[ , 2:5], method = "spearman"), 
            "clipboard",
            sep = "\t")

# ______________________________________________________________________________
# 6. Functional group composition ----

# here we'll show and summarize the implied composition by
  # - POO
  # - RRA

# It doesn't make much sense to use FOO because samples will be multi-counted
# BUT SSFOO sums to 1 so can't we use it?

# richness would be the same as POO in this context

# ______________________________________________________________________________
# 6a. Snow-off ----
# ______________________________________________________________________________

# total items
total.items.off <- length(unique(metrics.off$final.taxon))

# clean
df.for.plot.off <- bind_rows(
  
  # POO
  metrics.off |>
    
    filter(metric == "poo") |>
    
    group_by(func) |>
    summarize(comp = n() / total.items.off) |>
    
    mutate(metric = "poo"),
  
  # SSFOO
  metrics.off |>
    
    filter(metric == "wpoo") |>
    
    group_by(func) |>
    summarize(comp = sum(value)) |>
    
    mutate(metric = "ssfoo"),
  
  # RRA
  metrics.off |>
    
    filter(metric == "rra") |>
    
    group_by(func) |>
    summarize(comp = sum(value)) |>
    
    mutate(metric = "rra")
  
) |>
  
  # factor levels
  mutate(metric = factor(metric,
                         levels = c("poo", "ssfoo", "rra"),
                         labels = c("POO",
                                    "SSFOO",
                                    "RRA")))

# plot
plot.comp.off <- ggplot(data = df.for.plot.off) +
  
  theme_classic() +
  
  geom_bar(aes(x = metric,
               y = comp,
               fill = func),
           color = "white",
           stat = "identity") +
  
  scale_fill_manual(values = cat1.colors) +
  
  coord_cartesian(expand = F) +
  
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 360 - 45,
                                   vjust = 1,
                                   hjust = 0),
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     labels = c(0, 25, 50, 75, 100)) +
  
  ylab("Percent composition") +
  ggtitle("a")

plot.comp.off
  
# ______________________________________________________________________________
# 6b. Snow-on ----
# ______________________________________________________________________________

# total items
total.items.on <- length(unique(metrics.on$final.taxon))

# clean
df.for.plot.on <- bind_rows(
  
  # POO
  metrics.on |>
    
    filter(metric == "poo") |>
    
    group_by(func) |>
    summarize(comp = n() / total.items.on) |>
    
    mutate(metric = "poo"),
  
  # SSFOO
  metrics.on |>
    
    filter(metric == "wpoo") |>
    
    group_by(func) |>
    summarize(comp = sum(value)) |>
    
    mutate(metric = "ssfoo"),
  
  # RRA
  metrics.on |>
    
    filter(metric == "rra") |>
    
    group_by(func) |>
    summarize(comp = sum(value)) |>
    
    mutate(metric = "rra")
  
) |>
  
  # factor levels
  mutate(metric = factor(metric,
                         levels = c("poo", "ssfoo", "rra"),
                         labels = c("POO",
                                    "SSFOO",
                                    "RRA")))

# plot
plot.comp.on <- ggplot(data = df.for.plot.on) +
  
  theme_classic() +
  
  geom_bar(aes(x = metric,
               y = comp,
               fill = func),
           color = "white",
           stat = "identity") +
  
  scale_fill_manual(values = cat1.colors) +
  
  coord_cartesian(expand = F) +
  
  theme(axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 360 - 45,
                                   vjust = 1,
                                   hjust = 0),
        axis.title = element_blank(),
        legend.position = "none") +
  
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     labels = c(0, 25, 50, 75, 100),
                     position = "right") +
  
  ylab("Percent composition") +
  ggtitle("b")

plot.comp.on

# ______________________________________________________________________________
# 6c. Plot together ----
# ______________________________________________________________________________

plot_grid(plot.comp.off, plot.comp.on,
          rel_widths = c(1, 0.5))

# 574 x 380

# ______________________________________________________________________________
# 6d. Summaries ----
# ______________________________________________________________________________

df.comp.summ.off <- df.for.plot.off |>
  
  pivot_wider(values_from = comp,
              names_from = metric)

write.table(df.comp.summ.off, "clipboard", sep = "\t")

df.comp.summ.on <- df.for.plot.on |>
  
  pivot_wider(values_from = comp,
              names_from = metric)

write.table(df.comp.summ.on, "clipboard", sep = "\t")
