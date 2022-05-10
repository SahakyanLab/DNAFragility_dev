LoadKmertoneData <- function(
  breakpoint.experiment, 
  k, 
  kmer.table,
  kmer.ref.table,
  action = "ratio"
  ){
  
  # Parameter       Description
  # action          ratio, z-score, norm.control, case
  
  # load files
  files <- list.files(
    path = paste0("../data/kmertone/", breakpoint.experiment, "/"),
    pattern = paste0("score_", k)
  )
  
  data.set <- fread(
    file = paste0("../data/kmertone/", breakpoint.experiment, "/", files),
    sep = ",", header = TRUE, showProgress = FALSE,
    select = c("case", "control", "z")
  )
  
  # add kmer column
  data.set[, `:=` (kmer = kmer.table)]
  
  # only kep first occurring kmer in lexicological order 
  data.set <- data.set[match(kmer.ref.table$kmer, data.set$kmer)]
  
  # obtain weighting factors for use in breakpoint model
  if(action == "ratio"){
    norm.case <- data.set$case/sum(data.set$case, na.rm = TRUE)
    norm.control <- data.set$control/sum(data.set$control, na.rm = TRUE)
    ratio <- norm.case/norm.control 
    return(ratio)
  } else if(action == "z-score"){
    return(data.set$z)
  } else if(action == "case"){
    return(data.set$case)
  } else if(action == "norm.case"){
    norm.case <- data.set$case/sum(data.set$case, na.rm = TRUE)
    return(norm.case)
  } else if(action == "norm.control"){
    norm.control <- data.set$control/sum(data.set$control, na.rm = TRUE)
    return(norm.control)
  } 
}