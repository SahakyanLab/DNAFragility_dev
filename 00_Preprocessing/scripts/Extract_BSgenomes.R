# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

# to install this package
packages <- c(
    "BSgenome.Hsapiens.UCSC.hg18",
    "BSgenome.Hsapiens.UCSC.hg19",
    "BSgenome.Hsapiens.UCSC.hg38",
    "BSgenome.Hsapiens.1000genomes.hs37d5"
)

dir.to.create <- c(
  "hg18", "hg19", "hg38", "hs37d5", "1000_Genomes_exp"
)

for(package in packages){
    if(!package %in% rownames(installed.packages())){
        BiocManager::install(package)
    }
}

suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))

for(i in seq_along(dir.to.create)){
  if("1000_Genomes_exp" == dir.to.create[i]){
    system("wget https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz")
    system("mv human_g1k_v37.fasta.gz ../../Raw_data/ref/1000_Genomes_exp/")
    genome <- readDNAStringSet(filepath = "../../Raw_data/ref/1000_Genomes_exp/human_g1k_v37.fasta.gz")
  } else {
    suppressPackageStartupMessages(suppressWarnings(
      library(packages[i], character.only = TRUE)
    ))
  }

  # extract all pairs of chromosomes
  if(i == 1){
    genome <- BSgenome.Hsapiens.UCSC.hg18
  } else if(i == 2){
    genome <- BSgenome.Hsapiens.UCSC.hg19
  } else if(i == 3){
    genome <- BSgenome.Hsapiens.UCSC.hg38
  } else if(i == 4){
    genome <- BSgenome.Hsapiens.1000genomes.hs37d5
  }

  if(!dir.exists(paste0("../../Raw_data/ref/", dir.to.create[i], "/"))){
      dir.create(paste0("../../Raw_data/ref/", dir.to.create[i], "/"))
  }

  check.files.exist <- list.files(
    path = paste0("../../Raw_data/ref/", dir.to.create[i]),
    pattern = "*.fasta.gz"
  )

  if(length(check.files.exist) < 22){
    out <- pblapply(1:22, function(x){
      if("1000_Genomes_exp" == dir.to.create[i]){
        chr <- genome[i]
      } else {
        chr <- DNAStringSet(genome[[x]])
      }

      # save as fasta file
      writeXStringSet(
        x = chr,
        filepath = paste0("../../Raw_data/ref/", dir.to.create[i], 
                          "/chr", x, ".fasta.gz"), 
        format = "fasta",
        compress = TRUE
      )
    })
  }
  
  if("1000_Genomes_exp" == dir.to.create[i]){
    system("rm ../../Raw_data/ref/1000_Genomes_exp/human_g1k_v37.fasta.gz")
  }
}