# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
chromosome <- as.character(args[3])
ref.seq <- as.character(args[4])
k <- as.integer(args[5])
cores <- as.integer(args[6])
control <- as.logical(as.character(args[7]))

# my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/01_RMSD/scripts/"
# breakpoint.experiment="01-Nebulization/Pilot/1000_Genomes_exp_1"
# chromosome=1
# ref.seq="1000_Genomes_Pilot"
# k=4
# cores=1
# control=FALSE
setwd(my.path)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(pbapply))
if(length(args) > 0) pbo = pboptions(type="txt")

source("../lib/CalcKmerFreq.R")
source("../lib/CalcRMSD.R")
source("../../00_ReadsAlignment/lib/LoadBreakpoints.R")
# data.table::setDTthreads(threads = 1)

# import all breakpoint data
df <- LoadBreakpoints(
  experiment.folder = breakpoint.experiment, 
  chromosome = chromosome
)

if("lev.dist" %in% colnames(df)){
  df[, lev.dist := NULL]
}

# load reference sequence
ref.seq <- readDNAStringSet(filepath = paste0("../../Raw_data/ref/", 
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

# freq <- pbsapply(-50:51, function(x){
#   CalcKmerFreq(ind = x, k = k)$freq
# }, cl = cores)

rmsd.values <- pbsapply(1:(dim(freq)[2]-1), function(x){
  CalcRMSD(freq[, x], freq[, (x+1)])
})

saveRDS(
  object = rmsd.values,
  file = paste0("../data/", breakpoint.experiment,
                ifelse(control, "/control_rmsd_kmer_", "/rmsd_kmer_"),
                k, ".Rdata")
)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
limits <- length(rmsd.values)/2-1

p = as_tibble(rmsd.values) %>% 
  mutate(x = -limits:(length(rmsd.values)-limits-1)) %>% 
  ggplot(aes(x = x, y = value)) + 
  geom_point(size = 1) + 
  geom_vline(xintercept = -1.5:1.5, alpha = 0.3)

ggsave(filename = "Rplots.pdf", plot = p)