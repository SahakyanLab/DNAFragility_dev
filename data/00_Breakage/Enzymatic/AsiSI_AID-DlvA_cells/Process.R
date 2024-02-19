suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))

df <- fread("./AsiSI_restriction_sites.single.bed")

# AsiSI recognises GCGATCGC sites
df[, `:=`(V4 = NULL, V5 = NULL, V7 = NULL)]
setnames(df, c("seqnames", "start", "end", "strand"))
df <- plyranges::as_granges(df) %>% 
    dplyr::mutate(start = start + 1, end = end)

# extract sequences and ensure it's aligned correctly
hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
sequences <- getSeq(hg19, df)

# measure fragility at the middle of the sequence
df <- as_tibble(df) %>% 
    dplyr::mutate(start = start + 4, end = start, width = 1) %>% 
    dplyr::rename_with(~paste0("V", 1:5)) %>% 
    as.data.table()

fwrite(df, "output.bed")
unlink("./AsiSI_restriction_sites.single.bed")