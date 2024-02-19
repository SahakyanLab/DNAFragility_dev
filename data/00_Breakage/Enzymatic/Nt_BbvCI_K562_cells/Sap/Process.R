# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

sourceCpp("../../../../../02_Alignment/lib/edlibFunction.cpp")

# concatenate all files
files <- list.files(path = ".", pattern = ".bed", full.names = TRUE)
df <- lapply(files, fread)
df <- rbindlist(df) %>% as_tibble()
df <- df %>% 
    dplyr::filter(V1 %in% paste0("chr", 1:22)) %>%
    dplyr::arrange(V1, V2, V3) %>% 
    dplyr::select(V1, V2, V6) %>% 
    dplyr::mutate(V3 = 1, .after = V2) %>% 
    dplyr::rename_with(~c("seqnames", "start", "width", "strand")) %>% 
    dplyr::distinct() %>% 
    as.data.table()

# filter for breakages only in target region
# CC-TCAGC | GCTGA-GG
df <- plyranges::as_granges(df)
hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

dna.shift <- getSeq(
    hg19, 
    df %>% 
        dplyr::mutate(start = start+1, end = start) %>% 
        plyranges::anchor_3p() %>% 
        dplyr::mutate(width = 7)
)
query.seq <- paste0(dna.shift)

# forward
target.seq <- rep("GCTGAGG", length(query.seq))
lev.mat.fwd <- cbind(query.seq, target.seq)
lev.mat.fwd <- LevenshteinLoop(lev.mat.fwd)
lev.mat <- cbind(rev = lev.mat.fwd)
lev.mat <- apply(lev.mat, 1, min)

# discard rest
ind.of.int <- which(lev.mat <= 0)

# Nick bot / Nick top
# GC|TGAGG / CC|TCAGC
# CG|ACTCC / GG|AGTCG
breakpos <- getSeq(
    hg19,
    df[ind.of.int] %>%
        dplyr::filter(strand == "-") %>% 
        head(10) %>% 
        as_tibble() %>% 
        dplyr::select(-width) %>% 
        dplyr::mutate(
            start = ifelse(strand == "+", start-4, end+7),
            end = ifelse(strand == "+", end+1, end+7)
        ) %>% 
        plyranges::as_granges()
)

sequences <- getSeq(
    hg19, 
    df[ind.of.int] %>% 
        dplyr::mutate(start = start+1, end = start) %>% 
        plyranges::anchor_3p() %>% 
        dplyr::mutate(width = 7)
)


df_bp <- as_tibble(df[ind.of.int]) %>% 
    dplyr::mutate(start = ifelse(strand == "+", start-4, end+7)) %>%  
    dplyr::select(-end) %>%
    as.data.table()

dir.create(path = "./breakpoint_positions/", showWarnings = FALSE, recursive = TRUE)
dir.create(path = "./kmertone/", showWarnings = FALSE, recursive = TRUE)
for(chr in 1:22){
    temp <- df_bp[seqnames == paste0("chr", chr)]
    
    # breakpoints
    fwrite(
        data.table(start.pos = temp$start),
        paste0("./breakpoint_positions/chr", chr, ".csv")
    )

    # kmertone
    fwrite(
        data.table(
            chromosome = temp$seqnames,
            start.pos = temp$start,
            strand = temp$strand
        ),
        paste0("./kmertone/chr", chr, ".csv")
    )
}