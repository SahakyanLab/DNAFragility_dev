# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
ind <- as.numeric(args[2])
control <- as.logical(args[3])
seed <- as.numeric(args[4])
k <- as.numeric(args[5])

# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(cluster)))
suppressPackageStartupMessages(suppressWarnings(library(doParallel)))
suppressPackageStartupMessages(suppressWarnings(library(foreach)))
suppressPackageStartupMessages(suppressWarnings(library(mclust)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults

my.path="/Users/paddy/Documents/DPhil/03_Breakpoints_v2/03_FitCurves/scripts/"
setwd(my.path)
source("../lib/SequenceEffect.R")
source("../lib/FitCurves.R")

# calculate RMSD 
seq_effect <- SequenceEffect$new(
    chr = 1,
    which_exp_ind = ind,
    control = FALSE,
    cores = cores
)
seq_effect$calc_seq_effect(
    k = k, 
    rmsd.range = c(-301, 301)
)

# fit Gaussian curves
fit_curves <- FitCurves$new(
    chr = 1, 
    which_exp_ind = ind
)
fit_curves$generate_rmsd_plots()
