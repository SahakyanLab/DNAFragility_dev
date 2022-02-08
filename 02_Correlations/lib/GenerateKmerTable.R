GenerateKmerTable <- function(k, ref = TRUE){
  # create reference table of kmers
  k.mers <- do.call(data.table::CJ, rep(list(c("A", "C", "G", "T")), k))
  
  k.mers <- k.mers[, do.call(paste0, .SD)]
  
  if(ref){
    rev.comp <- as.character(reverseComplement(DNAStringSet(k.mers)))

    k.mer.ref <- data.table(
      'kmer' = k.mers, 
      'rev.comp' = rev.comp
      )[, `:=`(cond = ifelse(seq(1:nrow(.SD)) < match(kmer, rev.comp), 
      TRUE, ifelse(kmer == rev.comp, TRUE, FALSE)))][cond == TRUE, .(kmer)]

    return(k.mer.ref)
  } else {
    return(k.mers)
  }
}