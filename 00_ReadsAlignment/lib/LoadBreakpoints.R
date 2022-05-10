LoadBreakpoints <- function(
  experiment.folder, 
  chromosome,
  select.col = ""
  ){
  fetch.file <- paste0(
    "../../Raw_data/", experiment.folder,
    "/breakpoint_positions/chr", 
    chromosome, ".csv"
  )
  
  # obtain breakpoints
  if(nzchar(select.col)) {
    columns = select.col
  } else {
    columns = NULL
  }

  df <- fread(
    file = fetch.file,
    select = columns,
    showProgress = FALSE
  )
  
  return(df)
}