setwd("/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/00_ReadsAlignment/scripts/")

# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(pbapply))
pbo = pboptions(type="txt")

# load 1000 Genomes reference sequence
ref.seq = readDNAStringSet(filepath = "../../data/ref/1000_Genomes_exp/human_g1k_v37.fasta")

# extract all pairs of chromosomes
out <- pblapply(1:22, function(i){
  # load individual chromosome
  chr <- ref.seq[i]
  
  # save as fasta file
  writeXStringSet(
    x = chr, 
    filepath = paste0("../../data/ref/1000_Genomes_exp/chr", i, ".fasta.gz"),   
    format = "fasta",
    compress = TRUE
  )
})
