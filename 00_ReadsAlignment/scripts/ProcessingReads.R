# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.numeric(args[1])
breakpoint.experiment <- as.character(args[2])
experiment.num <- as.character(args[3])
chromosome <- as.character(args[4])
ref.path <- as.character(args[5])
ind <- as.numeric(args[6])
fasta.lines <- as.numeric(args[7])
interval <- as.numeric(args[8])
BAM <- as.logical(as.character(args[9]))
alignment.strands <- as.character(args[10])

my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/00_ReadsAlignment/scripts/"
breakpoint.experiment="03-Ancient_DNA/Chagyrskaya_Neandertal"
experiment.num=1
chromosome=1
ref.path="Simons_exp"
ind=0
fasta.lines=1000
BAM=FALSE
alignment.strands="both"
setwd(my.path)

# load dependencies
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Rcpp))
source("../lib/AlignReads.R")
sourceCpp("../lib/edlibFunction.cpp")

# load reference sequence
if(ref.path == "1000_Genomes_Pilot"){
  ref.seq.original <- list.files(path = paste0("../../data/ref/", ref.path), pattern = ".fa.gz$")
  ref.seq.original <- fasta.index(paste0("../../data/ref/", ref.path, "/", ref.seq.original))
  ref.seq.original <- readDNAStringSet(ref.seq.original[chromosome, ])
} else {
  ref.seq.original <- readDNAStringSet(filepath = paste0("../../data/ref/", ref.path, 
                                                  "/chr", chromosome, ".fasta.gz"))
}

for(i in 1:ind){
  # measure execution time 
  start.time <- Sys.time()

  output <- AlignReads(
    chromosome.nr=chromosome,
    ind=i,
    fasta.lines=fasta.lines,
    reads.path=paste0(breakpoint.experiment, "_", experiment.num),
    interval=interval,
    BAM=BAM,
    alignment.strands=alignment.strands
  )

  # output %>%
  #   fwrite(file = paste0(
  #     "../data/reads/", breakpoint.experiment, "_", experiment.num,
  #     "/breakpoint_positions/chr", chromosome.nr, 
  #     ifelse(alignment.strands == "plus", "/plus_alignment_", 
  #     ifelse(alignment.strands == "minus", "/minus_alignment_", 
  #     "/alignment_"), 
  #     "file_", ind,".txt"))
  #   )

  # measure execution time
  end.time   <- Sys.time()
  time.diff  <- signif(end.time-start.time, 3)
  time.units <- attr(time.diff, "units")
  cat("Execution time for run ", i, "/", ind,": ", time.diff, " ",  time.units, "\n")
  rm(output)
}