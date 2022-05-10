# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
chromosome <- as.numeric(args[3])
ref.path <- as.character(args[4])
ind <- as.numeric(args[5])
fasta.lines <- as.numeric(args[6])
interval <- as.numeric(args[7])
BAM <- as.logical(args[8])
alignment.strands <- as.character(args[9])

my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/00_ReadsAlignment/scripts/"
breakpoint.experiment="01-Nebulization/1000_Genomes_exp_2"
ref.path="1000_Genomes_exp"

# breakpoint.experiment="00-Ultrasonication/Simons_exp_1"
# ref.path="hs37d5"

# chromosome=1
# ind=0
# fasta.lines=11308072
# interval=1000
# BAM=FALSE
# alignment.strands="plus"

setwd(my.path)

# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Rcpp))
source("../lib/AlignReads.R")
sourceCpp("../lib/edlibFunction.cpp")

# load reference sequence
if(!file.exists(path = paste0("../../Raw_data/ref/", ref.path, 
                              "/chr", chromosome, ".fasta.gz"))){
  if(ref.path == "1000_Genomes_Pilot"){
    ref.seq.original <- list.files(path = paste0("../../Raw_data/ref/", ref.path), pattern = ".fa.gz$")
    ref.seq.original <- fasta.index(paste0("../../Raw_data/ref/", ref.path, "/", ref.seq.original))
    ref.seq.original <- readDNAStringSet(ref.seq.original[chromosome, ])

    writeXStringSet(
      x = ref.seq.original,
      filepath = paste0("../../Raw_data/ref/", ref.path, 
                        "/chr", chromosome, ".fasta.gz"),
      format = "fasta",
      compress = TRUE
    )
  }
} else {
  ref.seq.original <- readDNAStringSet(filepath = paste0("../../Raw_data/ref/", ref.path, 
                                                  "/chr", chromosome, ".fasta.gz"))
}

for(i in 0:(ind-1)){
  start.time <- Sys.time()

  output <- AlignReads(
    chromosome.nr=chromosome,
    ind=i,
    fasta.lines=fasta.lines,
    breakpoint.experiment=breakpoint.experiment,
    interval=interval,
    BAM=BAM,
    alignment.strands=alignment.strands
  )

  output %>%
    fwrite(file = paste0(
      "../../Raw_data/", breakpoint.experiment,
      "/breakpoint_positions/chr", chromosome, 
      ifelse(alignment.strands == "plus", "/plus_alignment_", 
      ifelse(alignment.strands == "minus", "/minus_alignment_", 
      "/alignment_")), 
      "file_", i,".txt"))

  end.time   <- Sys.time()
  time.diff  <- signif(end.time-start.time, 3)
  time.units <- attr(time.diff, "units")
  cat(paste0("Execution time for run ", i, "/", ind-1, " in ", time.diff, " ", time.units, "\n"))
  rm(output)
}