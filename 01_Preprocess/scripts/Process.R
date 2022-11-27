# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
load.and.process <- as.logical(args[2])

# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))

setwd(my.path)
source("../lib/PreprocessFiles.R")

# begin preprocessing steps
preprocess <- PreprocessFiles$new()
if(load.and.process){
    preprocess$load_and_process_files()
} else {
    preprocess$get_genomes()
}