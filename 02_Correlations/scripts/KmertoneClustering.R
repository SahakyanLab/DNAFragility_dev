# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
k <- as.integer(args[3])
upper.limit <- as.integer(args[4])
action <- as.character(args[5])
setwd(my.path)

# load dependencies
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(pbapply))
source("../lib/GenerateKmerTable.R")
source("../lib/LoadKmertoneData.R")

# kmer reference
sample.kmer <- data.table(kmer = GenerateKmerTable(k, ref = FALSE))
kmer.ref <- GenerateKmerTable(k, ref = TRUE)

# correlation coefficients of k-meric enrichments
all.exp <- 1:upper.limit
all.cor.coef <- pbsapply(all.exp, function(i){
  # study one
  exp_1 <- LoadKmertoneData(
    breakpoint.experiment  = breakpoint.experiment,
    experiment.num = i,
    k = k,
    kmer.table = sample.kmer$kmer,
    kmer.ref.table = kmer.ref,
    action = action
  )

  # correlation of study one with all other experiments
  compare.exp <- sapply(all.exp, function(x){
    if(x == i){
      cor.coef <- 1
    } else {
      exp_2 <- LoadKmertoneData(
        breakpoint.experiment  = breakpoint.experiment,
        experiment.num = x,
        k = k,
        kmer.table = sample.kmer$kmer,
        kmer.ref.table = kmer.ref,
        action = action
      )

      # remove any infinite values
      mat <- matrix(c(exp_1, exp_2), ncol = 2)
      mat <- mat[!rowSums(!is.finite(mat)), ]
      
      # calculate pearson correlation coefficient
      cor.coef <- cor(x = mat[, 1], y = mat[, 2], method = "pearson")
    }
    
    return(cor.coef)
  })
  return(compare.exp)
})

colnames(all.cor.coef) <- paste0("exp_", all.exp)
rownames(all.cor.coef) <- paste0("exp_", all.exp)
title.extract <- str_extract(string = experiment.folder, pattern = "[0-9].+(?=/)")

write.csv(
  x = all.cor.coef, 
  file = paste0("../data/kmertone/", title.extract, 
                "/correlation_coefficient_matrix_", 
                action,"_kmer_", k, ".csv"),
  row.names = FALSE
)