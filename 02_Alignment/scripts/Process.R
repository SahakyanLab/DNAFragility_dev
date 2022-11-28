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

my.path="/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints_v2/02_Alignment/scripts"
setwd(my.path)
source("../lib/ProcessReads.R")
sourceCpp("../lib/edlibFunction.cpp")

# begin preprocessing steps
preprocess <- ProcessReads$new(
    interval = 1000000, 
    which_exp_ind = 48
)
preprocess$process_reads()

df <- fread("../../data/org_file.csv")
df <- df[Processed == "FALSE" & `DSB Map` == "TRUE"]
ind=i=48
df[i]
# private=self=NULL
# private$org_file=df
# self$interval = NULL
# self$which_exp_ind = 1
# chr=1
# id=0