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

my.path="/Users/paddy/Documents/DPhil/03_Breakpoints_v2/04_KmericAnalysis/scripts"
setwd(my.path)
source("../lib/AnalyseKmers.R")
source("../lib/Kmertone/kmertone.R")