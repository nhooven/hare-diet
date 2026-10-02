# function - bootstrap metric
# this will yield bootstrap SEs and CIs for any given metric by recursively calling the function
# 95% CI is probably fine

boot_metric <- function (.reads,
                         .iter = 1000,
                         .metric = c("foo", "poo", "wpoo", "rra")) {
  
  # list of samples to boot over
  sample.IDs <- unique(.reads$Final.sample.ID)
  
  # foo
  if (.metric == "foo") {
    
    foo <- calc_foo(.reads)
    
    # bootstrap by sample
    boot.foo <- data.frame()
    
    #tic()
    
    for (i in 1:.iter) {
      
      samples.boot <- sample.IDs[sample(x = 1:length(sample.IDs),
                                        size = length(sample.IDs),
                                        replace = T)]
      
      # inner loop - "stack" samples and rename
      reads.samp <- data.frame()
      
      for (j in 1:length(samples.boot)) {
        
        reads.j <- .reads |> filter(Final.sample.ID == samples.boot[j]) |>
          
          # rename
          mutate(Final.sample.ID = j)
        
        reads.samp <- rbind(reads.samp, reads.j)
        
      } # j
      
      boot.foo <- rbind(boot.foo, calc_foo(reads.samp))
      
    } # i
    
    #toc() 
    # should take ~ 3 minutes for 1,000 boots
    
    # summarize
    boot.foo.sum <- boot.foo |>
      
      group_by(final.taxon) |>
      
      summarize(se = sd(pfoo),
                lci = quantile(pfoo, prob = 0.025),
                uci = quantile(pfoo, prob = 0.975))
    
    # join in
    df.out <- foo |> left_join(boot.foo.sum)
    
  }
  
  # poo
  if (.metric == "poo") {
    
    poo <- calc_poo(.reads)
    
    # bootstrap by sample
    boot.poo <- data.frame()
    
    #tic()
    
    for (i in 1:.iter) {
      
      samples.boot <- sample.IDs[sample(x = 1:length(sample.IDs),
                                        size = length(sample.IDs),
                                        replace = T)]
      
      # inner loop - "stack" samples and rename
      reads.samp <- data.frame()
      
      for (j in 1:length(samples.boot)) {
        
        reads.j <- .reads |> filter(Final.sample.ID == samples.boot[j]) |>
          
          # rename
          mutate(Final.sample.ID = j)
        
        reads.samp <- rbind(reads.samp, reads.j)
        
      } # j
      
      boot.poo <- rbind(boot.poo, calc_poo(reads.samp))
      
    } # i
    
    #toc() 
    # should take ~ 3 minutes for 1,000 boots
    
    # summarize
    boot.poo.sum <- boot.poo |>
      
      group_by(final.taxon) |>
      
      summarize(se = sd(poo),
                lci = quantile(poo, prob = 0.025),
                uci = quantile(poo, prob = 0.975))
    
    # join in
    df.out <- poo |> left_join(boot.poo.sum)
    
  }
  
  # wpoo
  if (.metric == "wpoo") {
    
    wpoo <- calc_wpoo(.reads)
    
    # bootstrap by sample
    boot.wpoo <- data.frame()
    
    #tic()
    
    for (i in 1:.iter) {
      
      samples.boot <- sample.IDs[sample(x = 1:length(sample.IDs),
                                        size = length(sample.IDs),
                                        replace = T)]
      
      # inner loop - "stack" samples and rename
      reads.samp <- data.frame()
      
      for (j in 1:length(samples.boot)) {
        
        reads.j <- .reads |> filter(Final.sample.ID == samples.boot[j]) |>
          
          # rename
          mutate(Final.sample.ID = j)
        
        reads.samp <- rbind(reads.samp, reads.j)
        
      } # j
      
      boot.wpoo <- rbind(boot.wpoo, calc_wpoo(reads.samp))
      
    } # i
    
    #toc() 
    # should take ~ 3 minutes for 1,000 boots
    
    # summarize
    boot.wpoo.sum <- boot.wpoo |>
      
      group_by(final.taxon) |>
      
      summarize(se = sd(wpoo),
                lci = quantile(wpoo, prob = 0.025),
                uci = quantile(wpoo, prob = 0.975))
    
    # join in
    df.out <- wpoo |> left_join(boot.wpoo.sum)
    
  }
  
  # rra
  if (.metric == "rra") {
    
    rra <- calc_rra(.reads)
    
    # bootstrap by sample
    boot.rra <- data.frame()
    
    #tic()
    
    for (i in 1:.iter) {
      
      samples.boot <- sample.IDs[sample(x = 1:length(sample.IDs),
                                        size = length(sample.IDs),
                                        replace = T)]
      
      # inner loop - "stack" samples and rename
      reads.samp <- data.frame()
      
      for (j in 1:length(samples.boot)) {
        
        reads.j <- .reads |> filter(Final.sample.ID == samples.boot[j]) |>
          
          # rename
          mutate(Final.sample.ID = j)
        
        reads.samp <- rbind(reads.samp, reads.j)
        
      } # j
      
      boot.rra <- rbind(boot.rra, calc_rra(reads.samp))
      
    } # i
    
    #toc() 
    # should take ~ 3 minutes for 1,000 boots
    
    # summarize
    boot.rra.sum <- boot.rra |>
      
      group_by(final.taxon) |>
      
      summarize(se = sd(rra),
                lci = quantile(rra, prob = 0.025),
                uci = quantile(rra, prob = 0.975))
    
    # join in
    df.out <- rra |> left_join(boot.rra.sum)
    
  }
  
  # return
  return(df.out)
  
}
