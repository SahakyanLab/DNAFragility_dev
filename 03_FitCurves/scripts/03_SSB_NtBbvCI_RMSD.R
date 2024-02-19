# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

#' 1. import breakpoint data
df <- fread("../../data/00_Breakage/Ultrasonication/Simons_exp_10/breakpoint_positions/chr17.csv")
df[, `:=`(seqnames = 17, strand = "+")]

ref <- Biostrings::readDNAStringSet(
    filepath = paste0("../../data/ref/", "hs37d5", 
                      "/chr", 17, ".fasta.gz")
)

hs37d5 <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
refseq.table <- as.data.frame(hs37d5@seqinfo)
refseq.table <- refseq.table[grepl(
    pattern = "^chr([1-9]|1[0-9]|2[0-2])$", 
    x = rownames(refseq.table)
),]
if(nrow(refseq.table) == 0){
    refseq.table <- as.data.frame(hs37d5@seqinfo)
    refseq.table <- refseq.table[grepl(
        pattern = "^([1-9]|1[0-9]|2[0-2])$", 
        x = rownames(refseq.table)
    ),]
    rownames(refseq.table) <- rownames(refseq.table)
}
chrs_len <- refseq.table$seqlengths         
chrs_len <- setNames(chrs_len, rownames(refseq.table))
rmsd.range <- -10:10

#' 2. make sure no position will be out-of-bounds
df[, `:=`(
    lower_end = start.pos-(max(abs(rmsd.range))+ceiling(kmer/2)),
    upper_end = start.pos+(max(abs(rmsd.range))+ceiling(kmer/2))
)]
df[, to_keep := ifelse(
    lower_end > 0 & upper_end < chrs_len["17"],
    TRUE, FALSE
)]
df <- df[to_keep == TRUE]
df[, `:=`(lower_end = NULL, upper_end = NULL, to_keep = NULL)]    

#' 3. generate reference table of 8-mers
k <- 8
k.mers <- do.call(data.table::CJ, 
                  rep(list(c("A", "C", "G", "T")), kmer))
kmer_list <- k.mers[, do.call(paste0, .SD)]
rev.comp <- as.character(
    Biostrings::reverseComplement(Biostrings::DNAStringSet(kmer_list))
)
kmer_ref <- data.table('fwd' = kmer_list, 'rev.comp' = rev.comp)
kmer_ref[, cond := ifelse(seq(1:nrow(.SD)) < match(fwd, rev.comp), 
TRUE, ifelse(fwd == rev.comp, TRUE, FALSE))]
kmer_ref <- kmer_ref[cond == TRUE, .(fwd, rev.comp)]

#' 4. calculate relative k-mer frequencies
calc_kmer_freq = function(ind, kmer){
    dfcopy <- copy(df)
    dfcopy[, `:=`(start.pos = start.pos + ind)]

    if((kmer %% 2) == 0){
        ending.pos <- ceiling((kmer-1)/2)
        starting.pos <- (kmer-1)-ending.pos
    } else {
        interval <- (kmer-1)/2
        starting.pos <- ending.pos <- interval
    }

    # obtain k-mers from alignment data
    dfcopy[,`:=`(start = start.pos - starting.pos, end = start.pos + ending.pos)]
    dfcopy[, start.pos := NULL]

    # extract k-meric counts and relative frequencies
    dfcopy[, `:=`(fwd = substring(text = ref, first = start, last = end))]

    if("-" %in% unique(dfcopy$strand)){
        dfcopy[strand == "-", fwd := paste(
            Biostrings::reverseComplement(Biostrings::DNAStringSet(fwd))
        )]
    }
    dfcopy <- dfcopy[, .(n = .N), by = .(fwd)][!stringr::str_detect(
        string = fwd, pattern = "N")]

    # account for strand symmetry
    fwd.ind <- match(kmer_ref$fwd, dfcopy$fwd)
    rev.comp.ind <- match(kmer_ref$rev.comp, dfcopy$fwd)

    # update data frame with dyad frequency count
    kmer_ref[, `:=`(freq = dfcopy$n[fwd.ind] + dfcopy$n[rev.comp.ind])]
    res <- kmer_ref$freq
    kmer_ref[, freq := NULL]

    # return data frame of k-mers with associated breakpoint frequencies
    return(res)
}

ind=0; k=k

freq <- pbapply::pblapply(rmsd.range, function(x){
  freq.vals <- CalcKmerFreq(ind = x, k = k)$freq
  norm.vals <- freq.vals/sum(freq.vals, na.rm = TRUE)
  return(list(freq.vals, norm.vals))
})

#' 5. calculate RMSD between adjacent positions
freq.vals <- sapply(freq, `[[`, 1)
norm.vals <- sapply(freq, `[[`, 2)

CalcRMSD <- function(a, b){
  return(sqrt(mean((a-b)^2, na.rm = TRUE)))
}

rmsd.values <- sapply(1:(dim(norm.vals)[2]-1), function(x){
  CalcRMSD(norm.vals[, x], norm.vals[, (x+1)])
})

#' 6. plot RMSD results
limits <- length(rmsd.values)/2-1

as_tibble(rmsd.values) %>% 
  mutate(x = -limits:(length(rmsd.values)-limits-1)) %>% 
  ggplot(aes(x = x, y = value)) + 
  geom_point(size = 1) + 
  geom_line()
