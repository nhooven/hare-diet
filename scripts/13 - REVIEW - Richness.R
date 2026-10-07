# PROJECT: Diet
# SCRIPT: 13 - REVIEW - Richness
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 06 Oct 2026
# COMPLETED: 07 Oct 2026
# LAST MODIFIED: 07 Oct 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(iNEXT)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# taxa information (from script 07)
reads.taxa <- read.csv("data_cleaned/reads_cleaned.csv")

# other studies
data.richness <- read.csv("data_review/data_richness.csv")

# ______________________________________________________________________________
# 3. Clean data ----
# ______________________________________________________________________________
# 3a. Case study ----
# ______________________________________________________________________________

# ROWS: taxa
# COLUMNS: samples

# ______________________________________________________________________________

# split by season
reads.off <- reads.taxa |> filter(Season == "Off")
reads.on <- reads.taxa |> filter(Season == "On")


# function to pivot into presence/absence data
format_inext <- function (.reads) {
  
  # data.frame
  taxa.df <- .reads |>
  
    # keep only necessary columns (sample, Final.taxon)
    dplyr::select(Final.sample.ID, final.taxon) |>
    
    # add a presence column
    mutate(presence = 1) |>
    
    # cast wide
    pivot_wider(names_from = Final.sample.ID,
                values_from = presence) |>
    
    # coerce everthing to a character
    mutate(across(everything(), as.character)) |>
    
    # replace "NULL" with 0, all else 1
    mutate(across(starts_with("S"), ~ifelse(. == "NULL",
                                            0,
                                            1))) |>
    
    # coerce to Boolean
    mutate(across(starts_with("S"), as.integer))
  
  # keep the first column as rownames, ensure the rest is an integer
  taxa.rows <- as.vector(taxa.df[ , 1])
  
  taxa.matrix <- as.matrix(taxa.df[ , -1])
  
  rownames(taxa.matrix) <- taxa.rows$final.taxon
  
  # return
  return(taxa.matrix)
  
}

# use function
off.matrix <- format_inext(reads.off)
on.matrix <- format_inext(reads.on)

# ______________________________________________________________________________
# 3b. Reviewed studies ----
# ______________________________________________________________________________

data.richness.1 <- data.richness |>
  
  group_by(Study) |>
  
  summarize(n.taxa = n()) |>
  
  ungroup() |>
  
  summarize(min = min(n.taxa),
            max = max(n.taxa))

# ______________________________________________________________________________
# 4. Rarefaction curves with iNEXT ----
# ______________________________________________________________________________

# define sequence of sample sizes
samp.seq <- 1:200

# off
off.iNEXT <- iNEXT(
  
  x = list("snow-off" = off.matrix),
  q = 0,
  datatype = "incidence_raw",
  size = samp.seq,
  se = T,
  conf = 0.95
  
)

# on
on.iNEXT <- iNEXT(
  
  x = list("snow-on" = on.matrix),
  q = 0,
  datatype = "incidence_raw",
  size = samp.seq,
  se = T,
  conf = 0.95
  
)

# prepare for plotting
forPlot <- bind_rows(fortify(off.iNEXT), fortify(on.iNEXT)) |>
  
  # inter/extrap
  mutate(prediction = ifelse(Method == "Extrapolation",
                             "extrapolation",
                             "interpolation")) |>
  
  mutate(prediction = factor(prediction,
                             levels = c("interpolation",
                                        "extrapolation")))

# ______________________________________________________________________________
# 5. Plot ----
# ______________________________________________________________________________

ggplot(data = forPlot) +
  
  theme_classic() +
  
  # min-max rect - snowshoe studies
  geom_rect(xmin = 0,
            xmax = 200,
            ymin = data.richness.1$min,
            ymax = data.richness.1$max,
            fill = "gray90") +
  
  geom_ribbon(aes(x = x,
                  y = y,
                  ymin = y.lwr,
                  ymax = y.upr,
                  fill = Assemblage),
              alpha = 0.25) +
  
  geom_line(aes(x = x,
                y = y,
                color = Assemblage,
                linetype = prediction),
            linewidth = 1) +
  
  # observed points
  geom_point(data = forPlot |> filter(Method == "Observed"),
             aes(x = x, 
                 y = y,
                 shape = Assemblage,
                 fill = Assemblage),
             size = 2) +
  
  # colors
  scale_color_manual(values = c("green4", "dodgerblue3")) +
  scale_fill_manual(values = c("green4", "dodgerblue3")) +
  scale_shape_manual(values = c(21, 23)) +
  
  # axis titles
  xlab("Number of samples") +
  ylab("Taxonomic richness") +
  
  # coords and scales
  coord_cartesian(expand = F) +
  scale_x_continuous(breaks = seq(25, 175, 25)) +
  
  # theme
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.9),
        legend.background = element_rect(fill = NA)) +
  
  # guides
  guides(linetype = "none") +
  
  # annotation
  annotate(geom = "text",
           x = 150,
           y = 12,
           size = 3,
           color = "gray30",
           label = "range from previous studies")

# 350 x 350

# NEXT
  # consider removing unknown taxa groups
  # split previous studies by season - worth being rigorous
