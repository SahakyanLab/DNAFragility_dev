# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])

# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
if(length(args) > 0) pbo <- pboptions(type = "txt")

my.path="/Users/paddy/Documents/DPhil/03_Breakpoints_v2/03_FitCurves/scripts"
setwd(my.path)
source("../lib/SequenceEffect.R")

# begin preprocessing steps
seq_effect <- SequenceEffect$new(
    chr = 1,
    which_exp_ind = 88,
    cores = 1,
    control = FALSE
)
seq_effect$calc_seq_effect(k = 4, rmsd.range = c(-21, 21))

# readRDS("../data/24-cfDNA/Bladder_cancer/rmsd_kmer_4.Rdata")

# df = fread("../../data/24-cfDNA/Bladder_cancer/breakpoint_positions/chr22.csv")
# df


# df <- fread("../../data/org_file.csv")
# df <- df[`DSB Map` == "TRUE"]
# private=self=NULL
# private$org_file=df
# self$chr=1
# self$which_exp_ind=88
# self$cores=1
# private$control=FALSE
# k=kmer=4
# ind=i=88
# rmsd.range=c(-10,10)

# x=0
# ind = x
# kmer = kmer