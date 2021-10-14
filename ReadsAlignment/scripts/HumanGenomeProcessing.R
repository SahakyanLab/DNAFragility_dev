setwd("/Volumes/Paddy_5TB/ProjectBoard_Patrick/03-Raw_Reads_Analysis/scripts/ReadsAlignment/")
suppressPackageStartupMessages(library(Biostrings))

# big files need more time to download
if(getOption('timeout') < 1000){
  options(timeout = 5000)
}

# download sequence
cat("Downloading human reference genome...")
if(!file.exists("../../data/GRCh38_latest_genomic.fna.gz")){
  download.file(url = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz", 
                destfile = "../../data/GRCh38_latest_genomic.fna.gz") 
}
cat("DONE!", "\n")

cat("Unpacking the toplevel combined fasta file...")
system(paste0("gunzip ../../data/GRCh38_latest_genomic.fna.gz"))
cat("DONE!", "\n")

# More memory efficient to request a small subset of sequences per iteration
fai <- fasta.index(paste0("../../data/GRCh38_latest_genomic.fna"))
fai <- readDNAStringSet(fai)

# extract base contents
cat("Extracting base contents...")
genome_length <- width(fai)
all.letters <- letterFrequency(fai, letters="ACGT", OR=0)
all.letters <- as.numeric(colSums(all.letters)/sum(genome_length))

# obtain other metadata
df <- data.frame(A.cont = all.letters[1],
                 C.cont = all.letters[2],
                 G.cont = all.letters[3],
                 T.cont = all.letters[4])

write.csv(df, file = "../../data/HumanGenomeATGC.csv", row.names = FALSE)
system("rm ../../data/GRCh38_latest_genomic.fna")
cat("DONE!", "\n")