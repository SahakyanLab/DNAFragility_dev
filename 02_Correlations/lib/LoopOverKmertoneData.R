LoopOverKmertoneData <- function(upper.limit, breakpoint.experiment, action){
  all.exp <- 1:upper.limit
  results <- pbsapply(all.exp, function(i){
    LoadKmertoneData(
      breakpoint.experiment = breakpoint.experiment,
      experiment = i,
      k = k,
      kmer.table = sample.kmer$kmer,
      kmer.ref.table = kmer.ref,
      action = action
    )
  })

  return(results)
}