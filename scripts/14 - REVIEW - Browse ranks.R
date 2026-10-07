# PROJECT: Diet
# SCRIPT: 14 - REVIEW - Browse ranks
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

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# case study
ranks.on <- read.csv("data_cleaned/summaries/ranks_on.csv")

# from reviewed studies
data.browse <- read.csv("data_review/data_browse.csv")

# ______________________________________________________________________________
# 3. Clean data ----
# ______________________________________________________________________________
# 3a. Case study -  Keep browse taxa only and re-rank ----
# ______________________________________________________________________________

ranks.browse <- ranks.on |>
  
  filter(final.taxon %in% c("Pinus contorta",
                            "Pinus albicaulis",
                            "Larix occidentalis",
                            "Alnus alnobetula",
                            "Shepherdia canadensis",
                            "Picea engelmannii",
                            "Abies lasiocarpa",
                            "Pseudotsuga menziesii",
                            "Salix scouleriana",
                            "Vaccinium scoparium",
                            "Juniperus communis",
                            "Ribes",
                            "Rhododendron neoglandosum",
                            "Paxistima myrsinites",
                            "Artemisia tridentata",
                            "Symphoricarpos",
                            "Populus",
                            "Rubus idaeus")) |>
  
  mutate(Genus = c("Pinus",
                   "Pinus",
                   "Larix",
                   "Alnus",
                   "Shepherdia",
                   "Picea",
                   "Abies",
                   "Pseudotsuga",
                   "Salix",
                   "Vaccinium",
                   "Juniperus",
                   "Ribes",
                   "Rhododendron",
                   "Paxistima",
                   "Artemisia",
                   "Symphoricarpos",
                   "Populus",
                   "Rubus")) |>
  
  # add rank columns
  mutate(
    
    # new ranking
    new.rank = 1:n(),
    
    # relative rank to mirror the other studies
    perc.rank = new.rank / max(new.rank)
    
    )

# extract Loomis genera
loomis.genera <- ranks.browse$Genus

# ______________________________________________________________________________
# 3b. Reviewed studies - "Relative" rank ----

# ranks are obviously contingent on how many taxa are included.
# we should re-calculate a "relative" rank (i.e., percent of total ranks)

# ______________________________________________________________________________

# carry study into study segment, if needed
data.browse <- data.browse |>
  
  # carry study into study segment, if needed
  mutate(Study_segment = ifelse(Study_segment == "",
                                Study,
                                Study_segment))

# bind in the max rank
data.browse.1 <- data.browse |>
  
  left_join(
    
    data.browse |>
      
      # group by study segment, determine max rank
      group_by(Study_segment) |>
      summarize(max.rank = max(importance_rank))
    
  ) |>
  
  # calculate percent of max rank
  mutate(perc.rank = importance_rank / max.rank)

# ______________________________________________________________________________
# 3c. Y-axis ordering - taxon-specific mean rank ----

# originally I included the case study here - really that's not what I want

# ______________________________________________________________________________

ordered.taxa.df <- data.browse.1 |>
  
  # keep only genus and perc.rank
  dplyr::select(Genus, perc.rank) |>
  
  # group by genus and calculate the mean rank
  group_by(Genus) |>
  summarize(mean.perc.rank = mean(perc.rank)) |>
  ungroup() |>
  
  # arrange
  arrange(desc(mean.perc.rank))

# vector
ordered.taxa <- ordered.taxa.df$Genus

# ______________________________________________________________________________
# 4. Base plot - all taxa ----
# ______________________________________________________________________________

ggplot() +
  
  theme_classic() +
  
  # lines for Loomis genera
  geom_hline(data = ranks.browse,
             aes(yintercept = Genus),
             alpha = 0.05,
             linewidth = 2,
             color = "dodgerblue4") +
  
  # reviewed studies - relative ranks
  geom_point(data = data.browse.1,
             aes(x = perc.rank,
                 y = Genus),
             shape = 21,
             color = "darkgray") +
  
  # case study - relative ranks
  geom_point(data = ranks.browse,
             aes(x = perc.rank,
                 y = Genus),
             shape = 23,
             size = 1.5,
             fill = "dodgerblue3") +
  
  # axes
  # correct axis ordering (this is way nicer!)
  scale_y_discrete(limits = ordered.taxa) +
  
  # reverse the axis
  scale_x_reverse(breaks = c(0.8, 0.2),
                  labels = c("low", "high")) +
  
  # titles
  xlab("Relative rank") +
  
  # theme
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic",
                                   size = 7),
        axis.ticks.x = element_blank())

# 429 x 529

# ______________________________________________________________________________
# 5. Focal plot - case study taxa ----
# ______________________________________________________________________________
# 5a. Filter data ----
# ______________________________________________________________________________

data.browse.2 <- data.browse.1 |> filter(Genus %in% loomis.genera)
ordered.taxa.2 <- ordered.taxa[ordered.taxa %in% loomis.genera]

# ______________________________________________________________________________
# 5b. Plot ----
# ______________________________________________________________________________

ggplot() +
  
  theme_classic() +
  
  # lines for Loomis genera
  geom_hline(data = ranks.browse,
             aes(yintercept = Genus),
             alpha = 0.05,
             linewidth = 2,
             color = "dodgerblue4") +
  
  # reviewed studies - relative ranks
  geom_point(data = data.browse.2,
             aes(x = perc.rank,
                 y = Genus),
             shape = 21,
             color = "darkgray") +
  
  # case study - relative ranks
  geom_point(data = ranks.browse,
             aes(x = perc.rank,
                 y = Genus),
             shape = 23,
             size = 1.5,
             fill = "dodgerblue3") +
  
  # axes
  # correct axis ordering (this is way nicer!)
  scale_y_discrete(limits = ordered.taxa.2) +
  
  # reverse the axis
  scale_x_reverse(breaks = c(0.8, 0.2),
                  labels = c("low", "high")) +
  
  # titles
  xlab("Relative rank") +
  
  # theme
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic",
                                   size = 7),
        axis.ticks.x = element_blank())

# 429 x 212

# ______________________________________________________________________________
# 6. Correlation ----
# ______________________________________________________________________________

# bind together
taxa.for.cor <- ordered.taxa.df |>
  
  left_join(
    
    ranks.browse |> dplyr::select(Genus, perc.rank),
    by = "Genus"
    
  ) |>
  
  # keep only complete cases
  drop_na()

# correlation
cor.test(taxa.for.cor$perc.rank, taxa.for.cor$mean.perc.rank)

# r = 0.624 [0.204, 0.850]

# plot
ggplot(data = taxa.for.cor) +
  
  theme_classic() +
  
  geom_point(aes(x = perc.rank,
                 y = mean.perc.rank),
             shape = 23,
             fill = "dodgerblue3") +
  
  geom_smooth(aes(x = perc.rank,
                  y = mean.perc.rank),
              method = "lm",
              color = "dodgerblue3",
              fill = "dodgerblue3",
              alpha = 0.15) +
  
  # labels
  ggrepel::geom_text_repel(aes(x = perc.rank,
                               y = mean.perc.rank,
                               label = Genus),
                           size = 2.5,
                           fontface = "italic") +
  
  scale_x_reverse(breaks = c(0.8, 0.2),
                  labels = c("low", "high")) +
  scale_y_reverse(breaks = c(0.8, 0.2),
                  labels = c("low", "high")) +
  
  coord_cartesian(xlim = c(0, 1.0),
                  ylim = c(0, 1.0)) +
  
  theme(axis.ticks = element_blank(),
        axis.text.y = element_text(angle = 90)) +
  
  xlab("Case study relative rank") +
  ylab("Reviewed studies relative rank")

# 350 x 350

# NEXT: 
  # compare all vs. just the "gold standard" selection studies?