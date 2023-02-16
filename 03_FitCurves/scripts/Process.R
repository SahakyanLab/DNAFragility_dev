# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])

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

my.path="/Users/paddy/Documents/DPhil/03_Breakpoints_v2/03_FitCurves/scripts"
setwd(my.path)
source("../lib/SequenceEffect.R")
source("../lib/FitCurves.R")
Rcpp::sourceCpp("../lib/SeqEffect.cpp")

# cores=ceiling(parallel::detectCores()/1.5)

# calculate RMSD 
seq_effect <- SequenceEffect$new(
    chr = 1,
    which_exp_ind = NULL,
    control = FALSE
)
seq_effect$calc_seq_effect(
    k = NULL, 
    rmsd.range = c(-501, 501)
)
seq_effect$quick_plot_check()
seq_effect$plots$kmer_4

# fit Gaussian curves
fit_curves <- FitCurves$new(
    chr = 1,
    which_exp_ind = NULL
)
fit_curves$generate_rmsd_plots()