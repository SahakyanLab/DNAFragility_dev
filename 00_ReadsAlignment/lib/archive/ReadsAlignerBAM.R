# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
chromosome.nr <- as.numeric(args[1])
ind <- as.numeric(args[2])
ref.path <- as.character(args[3])
reads.path <- as.character(args[4])
my.path <- as.character(args[5])
cores <- as.numeric(args[6])
fasta.lines <- as.numeric(args[7])
interval <- as.numeric(args[8])
alignment.strands <- as.character(args[9])

# chromosome.nr=1
# ind=0
# my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/00_ReadsAlignment/scripts/"
# ref.path="1000_Genomes_exp"
# reads.path="02-Sonication/1000_Genomes_exp_1"
# cores=1
# interval=10000
# alignment.strands="both"
setwd(paste0(my.path, "../lib/"))

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Rcpp))
sourceCpp(file = "edlibFunction.cpp")
pbo = pboptions(type="txt")

# measure execution time 
start.time <- Sys.time()

# load reference sequence
ref.seq.original <- readDNAStringSet(filepath = paste0("../../data/ref/", ref.path, 
                                                       "/chr", chromosome.nr, ".fasta.gz"))

if(ind == 0){
  reads <- fasta.index(paste0("../data/reads/", reads.path, "/chr", chromosome.nr, ".fasta.gz"),
                       skip = 0, 
                       nrec = interval)
} else {
  # load read sequences
  to.skip <- seq(from = interval+1, to = fasta.lines, by = interval+1) 
  reads <- fasta.index(paste0("../data/reads/", reads.path, "/chr", chromosome.nr, ".fasta.gz"),
                       skip = to.skip[ind], 
                       nrec = interval)
}
reads <- readDNAStringSet(reads)

# extract length of read
read.length <- width(reads)

# extract position on chromosome
read.names <- names(reads)
str.extract <- str_split(string = read.names, pattern = "_")

# strand information
strands <- sapply(str.extract, "[[", 2)

if(alignment.strands == "plus"){
  to.keep <- which(strands == "+")
  reads <- reads[to.keep]

  # bp start position
  read.names <- names(reads)
  read.length <- width(reads)
  str.extract <- str_split(string = read.names, pattern = "_")
  read.start.pos <- sapply(str.extract, "[[", 1)

  read.start.pos <- as.integer(read.start.pos)
  read.end.pos <- read.start.pos+read.length-1
  reads <- paste(reads)

  # get substrings of reference sequence
  mat <- matrix(c(read.start.pos, read.end.pos), ncol = 2)
  ref.seq <- str_sub(string = ref.seq.original, start = mat[,1], end = mat[,2])

  lev.dist <- pbsapply(1:length(reads), function(i){
    return(Levenshtein(seq1 = reads[i], seq2 = ref.seq[i]))
  })

  # concatenate results
  lev.dist <- t(lev.dist)
  lev.dist <- matrix(unlist(lev.dist[1,]), ncol = 1)

  # plus strand
  df <- data.frame(
    bp.start.pos = mat[,1]-1,
    lev.dist = lev.dist) %>%
    arrange(bp.start.pos) %>%
    group_by(bp.start.pos, lev.dist) %>%
    summarise(freq = n()) %>%
    filter(lev.dist == min(lev.dist)) %>%
    suppressMessages()

  df %>%
    fwrite(file = paste0("../data/reads/", reads.path, "/breakpoint_positions/chr",
                        chromosome.nr, "/plus_alignment_file_", ind,".txt"))

} else if(alignment.strands == "minus"){
  # reverse complement minus strands
  to.keep <- which(strands == "-")
  reads <- reverseComplement(reads[to.keep])

  # bp start position
  read.names <- names(reads)
  read.length <- width(reads)
  str.extract <- str_split(string = read.names, pattern = "_")
  read.start.pos <- sapply(str.extract, "[[", 1)

  read.start.pos <- as.integer(read.start.pos)
  read.end.pos <- read.start.pos+read.length-1
  reads <- paste(reads)

  # get substrings of reference sequence
  mat <- matrix(c(read.start.pos, read.end.pos), ncol = 2)
  ref.seq <- str_sub(string = ref.seq.original, start = mat[,1], end = mat[,2])

  lev.dist <- pbsapply(1:length(reads), function(i){
    return(Levenshtein(seq1 = reads[i], seq2 = ref.seq[i]))
  })

  # concatenate results
  lev.dist <- t(lev.dist)
  lev.dist <- matrix(unlist(lev.dist[1,]), ncol = 1)

  # plus strand
  df <- data.frame(
    bp.start.pos = mat[,2],
    lev.dist = lev.dist) %>%
    arrange(bp.start.pos) %>%
    group_by(bp.start.pos, lev.dist) %>%
    summarise(freq = n()) %>%
    filter(lev.dist == min(lev.dist)) %>%
    suppressMessages()

  df %>%
    fwrite(file = paste0("../data/reads/", reads.path, "/breakpoint_positions/chr",
                         chromosome.nr, "/minus_alignment_file_", ind,".txt"))
} else if(alignment.strands == "both"){
  to.revcomp <- which(strands == "-")
  reads <- c(reads[-to.revcomp], reverseComplement(reads[to.revcomp]))

  # bp start position
  read.names <- names(reads)
  read.length <- width(reads)
  str.extract <- str_split(string = read.names, pattern = "_")
  read.start.pos <- sapply(str.extract, "[[", 1)

  read.start.pos <- as.integer(read.start.pos)
  read.end.pos <- read.start.pos+read.length-1

  # get reverse complement of reads
  reads.revcomp <- reverseComplement(reads)
  reads.revcomp <- paste(reads.revcomp)
  reads <- paste(reads)

  # get substrings of reference sequence
  mat <- matrix(c(read.start.pos, read.end.pos), ncol = 2)
  ref.seq <- str_sub(string = ref.seq.original, start = mat[,1], end = mat[,2])

  lev.dist <- pbsapply(1:length(reads), function(i){
    return(list(Levenshtein(seq1 = reads[i], seq2 = ref.seq[i]),
                Levenshtein(seq1 = reads.revcomp[i], seq2 = ref.seq[i])))
  })

  # concatenate results
  lev.dist <- t(lev.dist)
  lev.dist <- matrix(c(unlist(lev.dist[,1]), unlist(lev.dist[,2])), ncol = 2)
  best.align <- max.col(-lev.dist)

  # keep read with lowest levenshtein distance
  which.plus.strand <- which(best.align == 1)
  plus.strand <- reads[which.plus.strand]

  which.minus.strand <- which(best.align == 2)
  minus.strand <- reads.revcomp[which.minus.strand]

  # concatenate best reads
  read <- vector(mode = "character", length = length(best.align))
  read[which.plus.strand] <- reads[which.plus.strand]
  read[which.minus.strand] <- reads.revcomp[which.minus.strand]

  # align first 2 nucleotides of read against reference sequence
  # plus strand
  kmer.read.plus <- str_sub(string = read[which.plus.strand], start = 1, end = 2)
  kmer.seq.plus <- str_sub(string = ref.seq[which.plus.strand], start = 1, end = 2)

  bp.start.pos.plus <- mat[which.plus.strand, 1]-1
  lev.dist.plus <- lev.dist[which.plus.strand, 1]

  kmer.seq.plus.match <- which(kmer.read.plus == kmer.seq.plus)
  bp.start.pos.plus <- bp.start.pos.plus[kmer.seq.plus.match]
  lev.dist.plus <- lev.dist.plus[kmer.seq.plus.match]

  # minus strand
  kmer.read.minus <- str_sub(string = read[which.minus.strand],
                            start = read.length[which.minus.strand]-1,
                            end = read.length[which.minus.strand])
  kmer.seq.minus <- str_sub(string = ref.seq[which.minus.strand],
                            start = read.length[which.minus.strand]-1,
                            end = read.length[which.minus.strand])

  bp.start.pos.minus <- mat[which.minus.strand, 2]
  lev.dist.minus <- lev.dist[which.minus.strand, 2]

  kmer.seq.minus.match <- which(kmer.read.minus == kmer.seq.minus)
  bp.start.pos.minus <- bp.start.pos.minus[kmer.seq.minus.match]
  lev.dist.minus <- lev.dist.minus[kmer.seq.minus.match]

  # combine results
  bp.start.pos <- c(bp.start.pos.plus, bp.start.pos.minus)
  lev.dist <- c(lev.dist.plus, lev.dist.minus)

  df <- data.frame(
    bp.start.pos = bp.start.pos,
    lev.dist = lev.dist) %>%
    arrange(bp.start.pos) %>%
    group_by(bp.start.pos, lev.dist) %>%
    summarise(freq = n()) %>%
    filter(lev.dist == min(lev.dist)) %>%
    suppressMessages()

  df %>%
    fwrite(file = paste0("../data/reads/", reads.path, "/breakpoint_positions/chr",
                        chromosome.nr, "/alignment_file_", ind,".txt"))
}

# measure execution time
end.time   <- Sys.time()
time.diff <- round(end.time-start.time, 3)
time.units <- attr(time.diff, "units")
cat(paste0("Execution time: ", time.diff, " ",  time.units, "\n"))