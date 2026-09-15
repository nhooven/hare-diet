# function - calculate proportion of occurrence
# this is the total number of sample-diet items which are a given taxon
calc_poo <- function (.reads) {
  
  # unique food items
  unique.items <- unique(.reads$final.taxon)
  n.items <- length(unique.items)
  
  # number of samples
  unique.samples <- unique(.reads$Final.sample.ID)
  n.samples <- length(unique.samples)
  
  # total sample-items (just the rows of .reads)
  n.sample.items <- nrow(.reads)
  
  # helper function - how many rows per item?
  calc_rows_per_item <- function (.item) {
    
    total.rows <- sum(.reads$final.taxon == .item)
    
  }
  
  # list and lapply
  unique.items.list <- as.list(unique.items)
  
  taxa.poo <- data.frame(
    
    final.taxon = unique.items,
    total.rows = do.call(rbind, lapply(unique.items.list, calc_rows_per_item))
    
  ) |>
    
    mutate(poo = total.rows / n.sample.items) |>
    
    # arrange by pfoo
    arrange(desc(poo))
  
  return(taxa.poo) 
  
}