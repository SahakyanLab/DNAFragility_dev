# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
cores <- as.numeric(args[3])
setwd(paste0(my.path, "../lib"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(pbapply))
pbo = pboptions(type="txt")
source("LoadBreakpoints.R")
# data.table::setDTthreads(threads = 1)

# ---------------
# obtain average levenshtein distance for all chromosomes
BreakpointsForKmertone <- pblapply(1:22, function(i){
  df <- LoadBreakpoints(
    experiment.folder = breakpoint.experiment, 
    chromosome = i,
    select.col = c("chromosome", "start.pos")
  )

  df[, chromosome := rep(paste0("chr", i), nrow(df))]
  setcolorder(df, c("chromosome", "start.pos"))
  
  fwrite(
    x = df, 
    row.names = FALSE, 
    file = paste0("../../Raw_data/", breakpoint.experiment, 
                  "/kmertone/chr", i, ".txt")
  )
}, cl = cores)