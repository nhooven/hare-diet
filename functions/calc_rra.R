# function - calculate relative read abundance
# this is the by-sample number of reads for each diet item, divided by the total reads

calc_rra <- function (.reads) {
  
  # total reads
  total.reads <- sum(.reads$reads)
  
  # number of samples
  unique.samples <- unique(.reads$Final.sample.ID)
  n.samples <- length(unique.samples)
  
  # calculate reads per sample
  calc_rel_reads <- function (.sample) {
    
    sample.reads <- .reads |> filter(Final.sample.ID == .sample) |>
      
      # calcualte relative reads
      mutate(rr = reads / sum(reads))
    
    return(sample.reads)
    
  }
  
  # list and apply
  unique.samples.list <- as.list(unique.samples)
  
  rel.reads <- do.call(rbind, lapply(unique.samples.list, calc_rel_reads))
  
  # sum weights by taxon, divide by total
  taxa.rra <- rel.reads |> group_by(final.taxon) |>
    
    summarize(all.rel.reads = sum(rr)) |>
    ungroup() |>
    mutate(rra = (all.rel.reads / n.samples)) |>
    dplyr::select(-all.rel.reads) |>
    arrange(desc(rra))
  
  return(taxa.rra)
  
}