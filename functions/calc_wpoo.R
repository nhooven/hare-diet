# function - calculate weighted proportion of occurrence
# this down-weights items in mixed-meal samples
# this is also called "split-sampled frequency of occurrence" - more intuitive

# Deagle et al 2019:
  # similar to POO, but rather than giving equal weight to all occurrences, 
  # this metric weights each occurrence according to the number of food items in the sample 
  # (e.g., if a sample contains five food items, each will be given weight 1/5)

calc_wpoo <- function (.reads) {
  
  # unique food items
  unique.items <- unique(.reads$final.taxon)
  n.items <- length(unique.items)
  
  # number of samples
  unique.samples <- unique(.reads$Final.sample.ID)
  n.samples <- length(unique.samples)
  
  # calculate down-weighted occurrence per sample
  calc_dwo <- function (.sample) {
    
    sample.reads <- .reads |> filter(Final.sample.ID == .sample) |>
      
      # add weight
      mutate(weight = 1 / n())
    
    return(sample.reads)
    
  }
  
  # list and apply
  unique.samples.list <- as.list(unique.samples)
  
  reads.with.weights <- do.call(rbind, lapply(unique.samples.list, calc_dwo))
  
  # sum weights by taxon, divide by total
  taxa.wpoo <- reads.with.weights |> group_by(final.taxon) |>
    
    summarize(total.weighted.occurrence = sum(weight)) |>
    ungroup() |>
    mutate(wpoo = total.weighted.occurrence / n.samples) |>
    dplyr::select(-total.weighted.occurrence) |>
    arrange(desc(wpoo))
  
  return(taxa.wpoo)
  
}