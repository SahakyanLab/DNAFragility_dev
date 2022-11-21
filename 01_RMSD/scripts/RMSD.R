# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
chromosome <- as.character(args[3])
ref.seq <- as.character(args[4])
k <- as.integer(args[5])
cores <- as.integer(args[6])
control <- as.logical(as.character(args[7]))
normalised <- as.logical(as.character(args[8]))

# my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/01_RMSD/scripts/"
# breakpoint.experiment="00-Ultrasonication/Simons_exp_1"
# chromosome=1
# ref.seq="hs37d5"
# k=4
# cores=1
# control=FALSE
setwd(my.path)

suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
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
setorder(df, start.pos)

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
# freq <- pblapply(-301:301, function(x){
#   freq.vals <- CalcKmerFreq(ind = x, k = k)$freq
#   norm.vals <- freq.vals/sum(freq.vals, na.rm = TRUE)
#   return(list(freq.vals, norm.vals))
# }, cl = cores)

freq <- pblapply(-50:51, function(x){
  freq.vals <- CalcKmerFreq(ind = x, k = k)$freq
  norm.vals <- freq.vals/sum(freq.vals, na.rm = TRUE)
  return(list(freq.vals, norm.vals))
}, cl = cores)

freq.vals <- sapply(freq, `[[`, 1)
norm.vals <- sapply(freq, `[[`, 2)

# save absolute frequency values for motif analysis
# saveRDS(
#   object = freq.vals,
#   file = paste0("../data/", breakpoint.experiment,
#                 ifelse(control, "/control_freq_rmsd_kmer_", "/freq_rmsd_kmer_"),
#                 k, ".Rdata")
# )

# save normalised frequency values
# saveRDS(
#   object = norm.vals,
#   file = paste0("../data/", breakpoint.experiment,
#                 ifelse(control, "/control_rmsd_kmer_", "/normalised_freq_rmsd_kmer_"),
#                 k, ".Rdata")
# )

rmsd.values <- pbsapply(1:(dim(norm.vals)[2]-1), function(x){
  CalcRMSD(norm.vals[, x], norm.vals[, (x+1)])
})

# saveRDS(
#   object = rmsd.values,
#   file = paste0("../data/", breakpoint.experiment,
#                 ifelse(control, "/control_rmsd_kmer_", "/rmsd_kmer_"),
#                 k, ".Rdata")
# )

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
limits <- length(rmsd.values)/2-1

p = as_tibble(rmsd.values) %>% 
  mutate(x = -limits:(length(rmsd.values)-limits-1)) %>% 
  ggplot(aes(x = x, y = value)) + 
  geom_point(size = 1) + 
  geom_line() + 
  geom_vline(xintercept = -1.5:1.5, alpha = 0.3)

ggsave(filename = "Rplots.pdf", plot = p)