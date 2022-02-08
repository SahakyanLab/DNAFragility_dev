# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
experiment.num <- as.character(args[3])
kmer <- as.integer(args[4])
setwd(paste0(my.path, "../lib/"))

# load dependencies
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))

# load files
files <- list.files(
  path = paste0("../data/kmertone/", breakpoint.experiment, "_", experiment.num), 
  pattern = paste0("score_", kmer)
)
files <- str_sort(files, numeric = TRUE)

# import files
data.set <- fread(
  file = paste0("../data/kmertone/", breakpoint.experiment, 
                "_", experiment.num, "/", files), 
  sep = ",",
  header = TRUE,
  showProgress = FALSE
)

# obtain weighting factors for use in breakpoint model
norm.case <- data.set$case/sum(data.set$case, na.rm = TRUE)
norm.control <- data.set$control/sum(data.set$control, na.rm = TRUE)
ratio <- norm.case/norm.control 

data.set <- data.table(
  kmer = data.set$kmer,
  weight.factor = ratio
)

fwrite(
  data.set,
  file = paste0("../data/weight_factor/", breakpoint.experiment, 
                "_", experiment.num, "/kmer_", kmer, ".csv"), 
  row.names = FALSE
)