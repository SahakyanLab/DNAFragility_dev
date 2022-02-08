# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
experiment.num <- as.character(args[3])
chromosome <- as.character(args[4])
ref.seq <- as.character(args[5])
k <- as.integer(args[6])
cores <- as.integer(args[7])
control <- as.logical(as.character(args[8]))

# my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/01_RMSD/scripts/"
# breakpoint.experiment="00-Ultrasonication/Simons_exp"
# experiment.num=1
# chromosome=1
# ref.seq="Simons_exp"
# k=2
# cores=2
# control=FALSE
setwd(my.path)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(pbapply))
# pbo = pboptions(type="txt")
source("../lib/CalcKmerFreq.R")
source("../lib/CalcRMSD.R")
source("../../00_ReadsAlignment/lib/LoadBreakpoints.R")
# data.table::setDTthreads(threads = 1)

# import all breakpoint data
df <- LoadBreakpoints(
  path.to.origin = "../../", 
  experiment.folder = breakpoint.experiment, 
  experiment = experiment.num, 
  chromosome = chromosome
)
df[, lev.dist := NULL]

# load reference sequence
ref.seq <- readDNAStringSet(filepath = paste0("../../data/ref/", 
                                              ref.seq, "/chr", chromosome, 
                                              ".fasta.gz"))

if(control){
  # create control breakpoints
  set.seed(seed = 1234)
  df <- data.table(
    start.pos = floor(runif(n = nrow(df), min = 1, max = width(ref.seq)))
    )[order(start.pos)]
}

# create reference table of kmers
k.mers <- do.call(data.table::CJ, rep(list(c("A", "C", "G", "T")), k))
k.mers <- k.mers[, do.call(paste0, .SD)]
rev.comp <- as.character(reverseComplement(DNAStringSet(k.mers)))
k.mer.ref <- data.table(
  'fwd' = k.mers, 
  'rev.comp' = rev.comp
  )[, `:=`(cond = ifelse(seq(1:nrow(.SD)) < match(fwd, rev.comp), 
  TRUE, ifelse(fwd == rev.comp, TRUE, FALSE)))][cond == TRUE, .(fwd, rev.comp)]

# RMSD calculations
freq <- pbsapply(-300:301, function(x){
  CalcKmerFreq(ind = x, k = k)$freq
}, cl = cores)

rmsd.values <- pbsapply(1:(dim(freq)[2]-1), function(x){
  CalcRMSD(freq[, x], freq[, (x+1)])
})

saveRDS(
  object = rmsd.values,
  file = paste0("../data/", breakpoint.experiment, "_", experiment.num,
                ifelse(control, "/control_rmsd_kmer_", "/rmsd_kmer_"),
                k, ".Rdata")
)
