FilterForExtremeKmers <- function(data, top = TRUE, group){
  if(top){
    kmers <- data[order(data[, get(group)], decreasing = FALSE)[1:5], "kmer"]
  } else {
    kmers <- data[order(data[, get(group)], decreasing = TRUE)[1:5], "kmer"]
  }
  
  kmers <- as.character(kmers[["kmer"]])
  return(kmers)
}