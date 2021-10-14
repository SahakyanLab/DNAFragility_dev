setwd("/Volumes/Paddy_5TB//ProjectBoard_Patrick/03-Raw_Reads_Analysis/scripts/ReadsAlignment/")

# to install this package
# BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

# This BSgenome data package was made from the following source data file:
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.1000genomes.hs37d5))

# progress bar
pb <- txtProgressBar(min = 1, max = 24, style = 3)

# extract all pairs of chromosomes
genome <- BSgenome.Hsapiens.1000genomes.hs37d5

for(i in 1:24){
  if(i == 23){
    chr <- DNAStringSet(genome[["X"]])
  } else if (i == 24) {
    chr <- DNAStringSet(genome[["Y"]])
  } else {
    chr <- DNAStringSet(genome[[i]])
  }
  
  # save as fasta file
  writeXStringSet(x = chr, 
                  filepath = paste0("../../data/ref/Simons_exp/chr", i, ".fasta"), 
                  format = "fasta")
  
  # update progress bar
  setTxtProgressBar(pb, i)
}