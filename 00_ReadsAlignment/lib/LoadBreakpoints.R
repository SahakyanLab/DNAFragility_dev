LoadBreakpoints <- function(
  path.to.origin = "./", 
  experiment.folder, 
  experiment, 
  chromosome,
  select.col = ""
  ){
  fetch.file <- paste0(path.to.origin, "00_ReadsAlignment/data/reads/",
                       experiment.folder, "_", experiment,
                       "/breakpoint_positions/chr", chromosome, ".csv")
  
  # obtain breakpoints
  df <- fread(
    file = fetch.file, 
    sep = ",", 
    header = TRUE, 
    select = ifelse(nzchar(select.col), select.col, NULL),
    showProgress = FALSE
  )
  
  return(df)
}