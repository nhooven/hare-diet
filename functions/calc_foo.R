# function - calculate (percent) frequency of occurrence
# this is the percent of all samples that include a given diet item

calc_foo <- function (.reads) {
  
  # unique food items
  unique.items <- unique(.reads$final.taxon)
  n.items <- length(unique.items)
  
  # number of samples
  unique.samples <- unique(.reads$Final.sample.ID)
  n.samples <- length(unique.samples)
  
  # helper function - check how many samples each taxon is in
  check_foo <- function (.taxon) {
    
    length(unique(.reads$Final.sample.ID[.reads$final.taxon == .taxon]))
    
  }
  
  # list and lapply
  unique.items.list <- as.list(unique.items)
  
  taxa.foo <- data.frame(
    
    final.taxon = unique.items,
    foo = do.call(rbind, lapply(unique.items.list, check_foo))
    
  ) |>
    
    mutate(pfoo = foo / n.samples) |>
    
    # arrange by pfoo
    arrange(desc(pfoo))
  
  return(taxa.foo)
  
}
