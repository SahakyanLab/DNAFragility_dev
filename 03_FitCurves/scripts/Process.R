args <- commandArgs(trailingOnly = TRUE)

# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(cluster)))
suppressPackageStartupMessages(suppressWarnings(library(mclust)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults

# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
k <- as.numeric(args[2])
setwd(my.path)

source("../lib/SequenceEffect.R")
source("../lib/FitCurves.R")
Rcpp::sourceCpp("../lib/SeqEffect.cpp")

# calculate RMSD 
seq_effect <- SequenceEffect$new(
    chr = 1:22,
    which_exp_ind = NULL,
    control = FALSE
)
seq_effect$calc_seq_effect(
    k = k, 
    rmsd.range = c(-501, 501)
)

# fit Gaussian curves
set.seed(seed)
fit_curves <- FitCurves$new(
    seed = 1234,
    chr = 1:22,
    which_exp_ind = NULL,
    return_vals = FALSE
)
fit_curves$generate_rmsd_plots(
    k = k, 
    per_chr = FALSE
)