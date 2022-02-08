# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
upper.limit <- as.integer(args[3])
setwd(paste0(my.path, "../lib"))

# load dependencies
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(gplots))

# import overlaps
nr.of.folders <- 1:upper.limit
path.of.experiment <- str_extract(string = breakpoint.experiment, pattern = "[0-9].+(?=/)")
files <- list.files(
  path = paste0("../data/overlap/", path.of.experiment), 
  pattern = "all_chromosomes.csv"
)

overlaps <- lapply(files, function(x){
  fread(
    file = paste0("../data/overlap/", path.of.experiment, "/", x), 
    sep = ",",
    header = TRUE,
    showProgress = FALSE
  )
})

overlaps <- lapply(overlaps, colMeans)

df <- lapply(overlaps, function(x){
  # which experiments are present
  df <- data.table(
    exp = names(x), 
    overlap = unname(x)
  )

  exp.in.df <- as.integer(str_extract(
      string = names(x), 
      pattern = "[0-9]+"
  ))
  
  # which experiment is missing
  exp.missing <- which(!nr.of.folders %in% exp.in.df)

  # append missing exp into df
  to.append <- data.table(
    exp = paste0("exp_", exp.missing),
    overlap = 1
  )

  df <- rbindlist(list(df, to.append))
  setorder(df, exp)
  return(df$overlap)
})

# concatenate results
df <- do.call(rbind, df)
colnames(df) <- paste0("exp_", nr.of.folders)
rownames(df) <- paste0("exp_", nr.of.folders)

# save overlap matrix
fwrite(
  as.data.frame(df),
  file = paste0("../data/overlap/", path.of.experiment, "/overlap_matrix.csv"),
  row.names = FALSE
)

pdf(paste0("../figures/", breakpoint.experiment, "_overlap_plot.pdf"))
heatmap.2(
  df,
  Rowv = FALSE,
  Colv = FALSE,
  dendrogram = "none",
  revC = FALSE,
  trace = "none",
  density.info = "none",
  cellnote = signif(df, 2),
  notecol = "black",
  notecex = 0.5,
  cexRow = 0.8,
  cexCol = 0.8,
  labRow = rownames(df),
  labCol = colnames(df)
)
dev.off()