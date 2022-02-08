setwd("/Volumes/Paddy_5TB/ProjectBoard_Patrick/03-Raw_Reads_Analysis/00_ReadsAlignment/scripts/")

# to install this package
# BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

# This BSgenome data package was made from the following source data file:
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.1000genomes.hs37d5))

# extract all pairs of chromosomes
genome <- BSgenome.Hsapiens.1000genomes.hs37d5

out <- pblapply(1:24, function(i){
  if(i == 23){
    chr <- DNAStringSet(genome[["X"]])
  } else if (i == 24) {
    chr <- DNAStringSet(genome[["Y"]])
  } else {
    chr <- DNAStringSet(genome[[i]])
  }
  
  # save as fasta file
  writeXStringSet(
    x = chr, 
    filepath = paste0("../../data/ref/Simons_exp/chr", i, ".fasta.gz"), 
    format = "fasta",
    compress = TRUE
  )
})