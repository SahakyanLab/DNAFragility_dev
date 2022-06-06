# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))

# obtain reference sequence
system("wget -r ftp://ftp.broadinstitute.org/pub/assemblies/mammals/elephant/loxAfr4/")

ref.folder <- "../../Raw_data/ref/savanna_elephant_loxAfr4/"
file.to.move <- "ftp.broadinstitute.org/pub/assemblies/mammals/elephant/loxAfr4/Chromosomes.v2.fasta.gz"

if(!dir.exists(ref.folder)){
    dir.create(ref.folder)
}

if(file.exists(file.to.move)){
    system(paste("mv", file.to.move, ref.folder))
    system("rm -r ftp.broadinstitute.org/")
    system(paste0("mv ", ref.folder, "Chromosomes.v2.fasta.gz ", ref.folder, "loxAfr4.fasta.gz"))
}

ref.files <- readDNAStringSet(filepath = paste0(ref.folder, "loxAfr4.fasta.gz"))

out <- pblapply(1:22, function(i){
  # save as fasta file
  writeXStringSet(
    x = ref.files[i], 
    filepath = paste0(ref.folder, "chr", i, ".fasta.gz"), 
    format = "fasta",
    compress = TRUE
  )
})

system(paste0("rm ", ref.folder, "loxAfr4.fasta.gz"))