# PROJECT: Diet
# SCRIPT: 17 - Metric diagram
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 08 Oct 2026
# COMPLETED: 08 Oct 2026
# LAST MODIFIED: 08 Oct 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 0. Purpose ----
# ______________________________________________________________________________

# demonstrate a toy visual example for each summary metric

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Overall parameters ----
# ______________________________________________________________________________

# total samples
n.samples <- 15

# total unique diet items
n.items <- 10

# ______________________________________________________________________________
# 3. Generate data ----

# each summary will be for one fictional plant species

# ______________________________________________________________________________

all.data <- expand.grid(
  
  sample = 1:n.samples,
  item = 1:n.items
  
)

# item 1 is our focus

# we'll start with richness per sample
rich.per.sample <- sample(1:n.items, size = n.samples, replace = T)

# which ones include item 1?
which.have1 <- rbinom(n = n.samples, size = 1, prob = 9/n.samples)

samp.with1 <- which(which.have1 == 1)

# richness minus item 1
rich.per.sample.minus1 <- rich.per.sample - which.have1

# list of which taxa
which.items <- list()

for (i in 1:n.samples) {
  
  # non-1 items
  which.items[[i]] <- sample(2:n.items, size = rich.per.sample.minus1[i])
  
  # add 1 if necessary
  if (i %in% samp.with1) {
    
    which.items[[i]] <- c(1, which.items[[i]])
    
  }
  
}

# subset original df
some.data <- data.frame()

for (i in 1:n.samples) {
  
  focal.data <- all.data |> filter(sample == i & item %in% which.items[[i]])
  
  some.data <- rbind(some.data, focal.data)
  
}

# add taxon
some.data <- some.data |>
  
  mutate(taxon = ifelse(some.data$item == 1, "focal", "other"))

# ______________________________________________________________________________
# 4. Plots ----
# ______________________________________________________________________________
# 4a. FOO ----
# ______________________________________________________________________________

data.foo <- data.frame(sample = 1:n.samples,
                       presence = which.have1)

plot.foo <- ggplot(data = data.foo) +
  
  theme_classic() +
  
  geom_col(aes(x = sample,
               y = presence),
           color = NA,
           fill = "dodgerblue3",
           linewidth = 0.2) +
  
  coord_cartesian(expand = F) +
  
  scale_x_continuous(breaks = 1:n.samples) +
  
  theme(axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 7)) +
  
  ylab("Occurrence")

plot.foo

# ______________________________________________________________________________
# 4b. POO ----
# ______________________________________________________________________________

plot.poo <- ggplot(data = some.data) +
  
  theme_classic() +
  
  geom_bar(aes(x = sample,
               fill = taxon,
               color = taxon,
               group = item),
           linewidth = 0.2,
           position = position_stack(reverse = T)) +
  
  coord_cartesian(expand = F) +
  
  scale_x_continuous(breaks = 1:n.samples) +
  
  theme(axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 7)) +
  
  scale_fill_manual(values = c("dodgerblue3", "lightgray")) +
  scale_color_manual(values = c(NA, "white")) +
  
  ylab("Count")

plot.poo

# ______________________________________________________________________________
# 4c. SSFOO ----
# ______________________________________________________________________________

plot.ssfoo <- ggplot(data = some.data) +
  
  theme_classic() +
  
  geom_bar(aes(x = sample,
               fill = taxon,
               color = taxon,
               group = item),
           linewidth = 0.2,
           position = position_fill(reverse = T)) +
  
  coord_cartesian(expand = F) +
  
  scale_x_continuous(breaks = 1:n.samples) +
  
  theme(axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 7)) +
  
  scale_fill_manual(values = c("dodgerblue3", "gray95")) +
  scale_color_manual(values = c(NA, "white")) +
  
  ylab("Weighted occurrence") 

plot.ssfoo

# ______________________________________________________________________________
# 4d. RRA ----
# ______________________________________________________________________________

# item-specific means
item.mean <- c(3000, runif(n.samples - 1, 150, 13000))

for (i in 1:nrow(some.data)) {
  
  some.data$reads[i] <- rpois(1, item.mean[some.data$item[i]] + rnorm(1, 0, 500))
  
}

hist(some.data$reads)

plot.rra <- ggplot(data = some.data) +
  
  theme_classic() +
  
  geom_col(aes(x = sample,
               y = reads,
               fill = taxon,
               color = taxon,
               group = item),
           linewidth = 0.2,
           position = position_fill(reverse = T)) +
  
  coord_cartesian(expand = F) +
  
  scale_x_continuous(breaks = 1:n.samples) +
  
  theme(axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 9)) +
  
  scale_fill_manual(values = c("dodgerblue3", "lightgray")) +
  scale_color_manual(values = c(NA, "white")) +
  
  xlab("Samples") +
  ylab("Relative reads")

plot.rra

# ______________________________________________________________________________
# 4e. Plot together ----
# ______________________________________________________________________________

cowplot::plot_grid(plot.foo, plot.poo, plot.ssfoo, plot.rra,
                   ncol = 1,
                   rel_heights = c(1, 1, 1, 1.2))

# 246 x 461

# ______________________________________________________________________________
# 5. Calculate metrics for demonstration ----
# ______________________________________________________________________________

# FOO
sum(which.have1) / n.samples  # 9 / 15 = 60%

# POO
length(which(some.data$item == 1)) / nrow(some.data)  # 9 / 65 = 13.8%

# SSFOO
data.ssfoo <- some.data |>
  
  group_by(sample) |>
  summarize(total.items = n()) |>
  ungroup() |>
  
  mutate(weight = 1 / total.items) |>
  
  right_join(some.data) |>
  
  filter(item == 1)

sum(data.ssfoo$weight) / n.samples # 2.16 / 15 = 0.144

# RRA
data.rra <- some.data |>
  
  group_by(sample) |>
  summarize(total.reads = sum(reads)) |>
  ungroup() |>
  
  right_join(some.data) |>
  
  mutate(rra = reads / total.reads) |>
  
  filter(item == 1)

rra.vec <- c(data.rra$rra, rep(0, 6))

sum(rra.vec) / n.samples   # 0.064
