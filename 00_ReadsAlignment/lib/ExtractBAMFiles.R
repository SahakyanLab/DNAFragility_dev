# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
setwd(paste0(my.path, "../lib"))

suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
pbo = pboptions(type="txt")

LoadChr <- function(seqname, bamFile, with.index){
  # reference names and lengths
  bamInfo <- seqinfo(bamFile) 
  
  # start and end region to extract
  afrom <- 1 
  ato <- seqlengths(bamInfo[as.character(seqname)]) 
  
  if(with.index){
    # set parameters for loading
    param <- ScanBamParam(
      what = c('seq', 'pos', 'strand'),
      which = GRanges(seqname, IRanges(afrom, ato)),
      flag = scanBamFlag(isUnmappedQuery = FALSE)
    )
    
    # load bam with region limits
    bam.content <- scanBam(file = bamFile, param=param)[[1L]] 
  } else {
    # load bam with region limits
    bam.content <- scanBam(file = bamFile)[[1L]] 
  }
  
  # extract sequence data
  sequences <- bam.content$seq
  names(sequences) <- paste0(bam.content$pos, "_", bam.content$strand)
  
  # save as fasta file
  writeXStringSet(
    x = sequences, 
    filepath = paste0("../../Raw_data/", breakpoint.experiment, 
                      "/chr", as.character(seqname), ".fasta.gz"), 
    format = "fasta", 
    compress = TRUE
  )
}

# locate BAM file
files <- list.files(
  path = paste0("../../Raw_data/", breakpoint.experiment),
  pattern = "*.bam"
)
files <- str_sort(files, numeric = TRUE)

if(length(files) > 2){
  files <- list.files(
    path = paste0("../../Raw_data/", breakpoint.experiment),
    pattern = "*.bam$"
  )
  files <- str_sort(files, numeric = TRUE)
  
  for(i in 1:length(files)){
    bamFile <- BamFile(paste0("../../Raw_data/", breakpoint.experiment, "/", files[i]))

    bam <- LoadChr(seqname = i, bamFile = bamFile, with.index = TRUE)
  }
} else {
  # check if index file exists
  check.ind.exist <- str_extract(string = files, pattern = "\\.([[:alnum:]]+)$")
  with.index <- ifelse(".bai" %in% check.ind.exist, TRUE, FALSE)

  bamFile <- BamFile(paste0("../../Raw_data/", breakpoint.experiment, "/", files[1]))

  seqnames <- 1:22
  bam <- pblapply(seqnames, LoadChr, bamFile, with.index)

  # Cleanup
  close(bamFile)
}