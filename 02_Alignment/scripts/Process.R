# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])

# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(Rsamtools)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))

my.path="/Users/paddy/Documents/DPhil/03_Breakpoints_v2/02_Alignment/scripts"
setwd(my.path)
source("../lib/ProcessReads.R")
sourceCpp("../lib/edlibFunction.cpp")

# begin preprocessing steps
preprocess <- ProcessReads$new(
    interval = 1000000, 
    all_exp = FALSE,
    which_exp_ind = 48
)
preprocess$process_reads()

df <- fread("../../data/org_file.csv")
df <- df[Processed == "FALSE" & `DSB Map` == "TRUE"]
ind=i=48
df[i]
private=self=NULL
private$org_file=df
self$interval = 1000000 
self$all_exp = FALSE
self$which_exp_ind = 48
chr=1
id=0

# breakpoint.experiment="24-cfDNA/Bladder_cancer"
# ref.path="hs37d5"
# chromosome=1
# ind=0
# fasta.lines=11308072
# interval=5000000
# BAM=FALSE
# BED=TRUE
# alignment.strands="both"

# reads <- fasta.index(
#     paste0("../../Raw_data/", breakpoint.experiment, 
#             "/chr", chromosome.nr, ".fasta.gz"),
#     skip = 0, 
#     nrec = interval
# )

# fasta.lines <- system(
#     command = "/usr/bin/zgrep -Ec '>' ../../data/24-cfDNA/Bladder_cancer/chr1.fasta.gz",
#     intern = TRUE
# )
# fasta.lines <- as.numeric(fasta.lines)
# ind <- floor(((fasta.lines+interval-1)/interval)-1)