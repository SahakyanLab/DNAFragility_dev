args <- commandArgs(trailingOnly = TRUE)

# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(gplots)))
suppressPackageStartupMessages(suppressWarnings(library(dendextend)))
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))
pbo <- pboptions(type = 'txt', char = '=')
data.table::setDTthreads(threads = 1) # prevents segmentation faults

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

k <- as.numeric(args[2])
action <- as.character(args[3])

source("../lib/Kmertone/kmertone.R")
source("../lib/AnalyseKmers.R")
Rcpp::sourceCpp("../lib/edlibFunction.cpp")

analyse_kmers <- AnalyseKmers$new(
    which_exp_ind = NULL,
    cores = 1,
    statistic = "mean"
)
analyse_kmers$get_kmer_enrichment(k = k)

analyse_kmers$get_kmer_tracts(
    k = k,
    action = action, 
    get_table = TRUE
)