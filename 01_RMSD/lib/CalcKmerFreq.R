CalcKmerFreq <- function(ind, k){
  dfcopy <- copy(df)
  dfcopy[, `:=`(start.pos = start.pos + ind)]
  
  if((k %% 2) == 0){
    ending.pos = ceiling((k-1)/2)
    starting.pos = (k-1)-ending.pos
  }

  # obtain k-mers from alignment data
  dfcopy[,`:=`(start = start.pos - starting.pos, end = start.pos + ending.pos)]
  dfcopy[, start.pos := NULL]

  # extract k-meric counts and relative frequencies
  dfcopy[, `:=`(fwd = str_sub(string = ref.seq, start = start, end = end))]

  if("-" %in% unique(dfcopy$strand)){
    dfcopy[strand == "-", fwd := paste(reverseComplement(DNAStringSet(fwd)))]
  }

  dfcopy <- dfcopy[, .(n = .N), by = .(fwd)][!str_detect(string = fwd, pattern = "N")]
  
  # account for strand symmetry
  fwd.ind <- match(k.mer.ref$fwd, dfcopy$fwd)
  rev.comp.ind <- match(k.mer.ref$rev.comp, dfcopy$fwd)
  
  # update data frame with dyad frequency count
  k.mer.ref[, `:=`(freq = dfcopy$n[fwd.ind] + dfcopy$n[rev.comp.ind])]
  k.mer.ref[, `:=`(freq = freq/sum(freq, na.rm = TRUE))]

  # return data frame of k-mers with associated breakpoint frequencies
  return(k.mer.ref[, .(freq)])
}