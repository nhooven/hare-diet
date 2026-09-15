# PROJECT: Diet
# SCRIPT: 10 - Visualize samples
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

# samples
samples.lookup <- read.csv("data_cleaned/samples_lookup.csv")

# cleaned reads
clean.reads <- read.csv("data_cleaned/reads_cleaned.csv")

# ______________________________________________________________________________
# 3. Clean ----
# ______________________________________________________________________________

samples.lookup.cut1 <- samples.lookup |> 
  filter(include.reasonable == "Y") |> 
  dplyr::select(-c(include.reasonable, include.stringent))

samples.lookup.cut2 <- samples.lookup |> 
  filter(include.stringent == "Y") |> 
  dplyr::select(-c(include.reasonable, include.stringent))

# keep only necessary reads
clean.reads.cut1 <- clean.reads |>
  filter(Final.sample.ID %in% samples.lookup.cut1$Final.sample.ID)

# seasonal splits
clean.reads.off <- clean.reads.cut1 |> filter(Season == "Off")
clean.reads.on <- clean.reads.cut1 |> filter(Season == "On") 

# ______________________________________________________________________________
# 4. Collection dates ----
# ______________________________________________________________________________
# 4a. Prep data ----
# ______________________________________________________________________________

samples.date <- samples.lookup.cut1 |>
  
  mutate(Collection.date = ymd(Collection.date),
         Collection.doy = yday(Collection.date))

# first of month lookup
month.lookup <- data.frame(
  
  date = paste0(c(paste0(0, 1:9), 10:12), "-01-2024")
  
) |>
  
  mutate(
    
    # doy
    doy = yday(mdy(date)),
    
    # label
    label = paste0(c("Jan", "Feb", "Mar", 
                     "Apr", "May", "Jun", 
                     "Jul", "Aug", "Sep", 
                     "Oct", "Nov", "Dec"),
                   " 01"))

# snow-on/off dates
snow.lookup <- data.frame(
  
  date = c("11-06-2022",
           "04-30-2023",
           "10-25-2023",
           "04-24-2024",
           "10-30-2024",
           "04-30-2025")
  
) |>
  
  mutate(doy = yday(mdy(date)))
    
# ______________________________________________________________________________
# 4b. Plot ----
# ______________________________________________________________________________

ggplot(samples.date) +
  
  theme_minimal() +
  
  # snow dates
  geom_vline(xintercept = snow.lookup$doy,
             linetype = "dashed") +
  
  geom_histogram(aes(x = Collection.doy,
                     fill = Season),
                 color = "black",
                 bins = 50) +
  
  # circular
  coord_radial(direction = -1,
               reverse = T,
               r.axis.inside = T) +
  
  # theme
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none") +
  
  # radial axis
  scale_y_continuous(breaks = c(10, 15, 20)) +
  
  # circumferential axis
  scale_x_continuous(breaks = month.lookup$doy,
                     labels = month.lookup$label,
                     expand = c(0, 0)) +
  
  # colors
  scale_fill_manual(values = c("green4", "dodgerblue2"))

# ______________________________________________________________________________
# 5. Sample richness ----
# ______________________________________________________________________________
# 5a. Prep data ----
# ______________________________________________________________________________

# function - rearranges cat2 functional groups
prep_rich <- function (.reads) {
  
  reads.out <- .reads |>
    
    group_by(Final.sample.ID, final.taxon) |>
    
    slice(1) |>
    
    ungroup() |>
    
    mutate(func = factor(cat2,
                         levels = c("lodgepole pine",
                                    "conifer",
                                    "woody broadleaf",
                                    "sub-shrub",
                                    "forb", 
                                    "graminoid",
                                    "unknown")))
  
  return(reads.out)
  
}

# use
clean.reads.off.1 <- prep_rich(clean.reads.off)
clean.reads.on.1 <- prep_rich(clean.reads.on)

# colors
cat2.colors <- c("darkblue",
                 "darkgreen", "darkorange3", "brown",
                 "green3", "lightgreen",
                 "gray")

# ______________________________________________________________________________
# 5b. Function ----
# ______________________________________________________________________________

plot_rich <- function(.reads) {
  
  ggplot(.reads) +
  
    theme_classic() +
    
    geom_bar(aes(x = Final.sample.ID,
                 fill = func),
             color = "white") +
    
    theme(legend.position = "right",
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    
    scale_y_continuous(expand = c(0, 0)) +
    
    scale_fill_manual(values = cat2.colors) +
    
    xlab("Samples") +
    ylab("Total diet items")
  
}

# ______________________________________________________________________________
# 5c. Plot ----
# ______________________________________________________________________________

plot_rich(clean.reads.off.1)
plot_rich(clean.reads.on.1)

# ______________________________________________________________________________
# 6. Sample RRA ----
# ______________________________________________________________________________
# 6a. Prep data ----
# ______________________________________________________________________________

# function - rearranges cat2 functional groups
prep_rra <- function (.reads) {
  
  reads.1 <- .reads |>
    
    # item reads
    group_by(Final.sample.ID, final.taxon) |>
    summarize(total.item.reads = sum(reads)) |>
    ungroup() |>
    
    # add in identifiers
    left_join(
      
      .reads |> dplyr::select(
        
        Final.sample.ID,
        final.taxon,
        cat1,
        cat2,
        Treatment,
        Season,
        Sex,
        Site,
        batch
        
      )
      
    ) |>
    
    mutate(func = factor(cat2,
                         levels = c("lodgepole pine",
                                    "conifer",
                                    "woody broadleaf",
                                    "sub-shrub",
                                    "forb", 
                                    "graminoid",
                                    "unknown"))) |>
    
    # drop duplicates
    distinct()
  
  # divide by total reads
  reads.out <- reads.1 |>
    
    group_by(Final.sample.ID) |>
    summarize(total.sample.reads = sum(total.item.reads)) |>
    ungroup() |>
    
    right_join(reads.1) |>
    
    # RRA
    mutate(rra = total.item.reads / total.sample.reads)
  
  return(reads.out)
  
}

# use
reads.off.rra <- prep_rra(clean.reads.off)
reads.on.rra <- prep_rra(clean.reads.on)

# ______________________________________________________________________________
# 6b. Function ----
# ______________________________________________________________________________

plot_rra <- function(.reads) {
  
  ggplot(.reads) +
    
    theme_classic() +
    
    geom_col(aes(x = Final.sample.ID,
                 y = rra,
                 fill = func),
             color = "white") +
    
    theme(legend.position = "right",
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    
    scale_y_continuous(expand = c(0, 0)) +
    
    scale_fill_manual(values = cat2.colors) +
    
    xlab("Samples") +
    ylab("Relative read abundance")
  
}

# ______________________________________________________________________________
# 5c. Plot ----

# 800 x 300

# ______________________________________________________________________________

plot_rra(reads.off.rra)
plot_rra(reads.on.rra)

# NEXT: 
  # thinner lines?
  # orde by collection date?
