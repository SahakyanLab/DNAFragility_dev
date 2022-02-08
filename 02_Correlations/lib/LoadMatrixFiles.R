LoadMatrixFiles <- function(action){
  files <- list.files(
    path = paste0("../data/kmertone/", title.extract, "/"), 
    pattern = paste0("correlation_coefficient_matrix_", action)
  )
  files <- str_sort(files, numeric = TRUE)

  data.sets <- sapply(files, function(x){
    # load matrices
    mat <- fread(
      file = paste0("../data/kmertone/", title.extract, "/", x),
      header = TRUE, 
      showProgress = FALSE
    )

    # flatten to numeric vector
    rownames(mat) <- colnames(mat)
    mat[lower.tri(mat, diag = TRUE)] <- NA
    mat <- as.numeric(na.omit(unlist(mat)))

    return(mat)
  })
  return(data.sets)
}