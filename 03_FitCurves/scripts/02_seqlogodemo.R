suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(ggseqlogo)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

#' 1. import breakage data for chromosome 1 of Nt.BbvCI and EcoRV enzymes
get_seqlogos <- function(f){
    df <- fread(f)
    setnames(df, "start")

    #' 2. expand into 8-mer
    df[, `:=`(seqnames = "chr1", width = 1, strand = "+")]
    setcolorder(df, c("seqnames", "start", "width", "strand"))
    df.granges <- plyranges::as_granges(df)

    df.expand <- plyranges::stretch(df.granges, 6L) 
    df.expand <- df.expand %>% mutate(end = end+1)

    #' 3. extract sequence from reference genome
    hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    df.expand.seq <- getSeq(hg19, df.expand)

    #' 4. get sequence logo
    df.expand.chars.fwd <- paste0(df.expand.seq)
    df.expand.chars.rc <- paste0(reverseComplement(df.expand.seq))
    df.expand.chars.both <- c(df.expand.chars.fwd, df.expand.chars.rc)

    seqlogo <- ggseqlogo(df.expand.chars.both, method = "bits", seq_type = 'dna', col_scheme = "nucleotide") + 
        theme(
            axis.line.x = element_line(colour = "black"),
            axis.ticks.x = element_line(colour = "black"),
            axis.line.y = element_line(colour = "black"),
            axis.ticks.y = element_line(colour = "black"),
            text = element_text(size = 15)
        ) + 
        coord_cartesian(ylim = c(0, 2)) %>% 
        suppressWarnings()

    return(seqlogo)
}

enz.nt <- paste0(
    "../../data/00_Breakage/Enzymatic/", 
    "Nt_BbvCI_K562_cells/NT/breakpoint_positions/",
    "chr1.csv"
)
get_seqlogos(f = enz.nt)

enz.ecorv <- paste0(
    "../../data/00_Breakage/Enzymatic/",
    "EcoRV_HeLa_cells/breakpoint_positions/", 
    "chr1.csv"
)
get_seqlogos(f = enz.ecorv)