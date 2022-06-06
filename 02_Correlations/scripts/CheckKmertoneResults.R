my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/02_Correlations/scripts/"
# breakpoint.experiment="10-DSBCapture/NHEK_DSBs"
breakpoint.experiment="04-Enzymatic/AID-DlvA_AsiSI_restriction_enzyme"
# breakpoint.experiment="19-CCseq/RPE-1_Top2-linked_DNA_DSBs/WG13-WG21_RPE1_WT_G1_UT"
breakpoint.experiment="25-RAFT/DSBs_in_HEK293T"
chromosome=1
k=4
setwd(my.path)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(library(ggplot2))

source("../../00_ReadsAlignment/lib/LoadBreakpoints.R")
source("../lib/GenerateKmerTable.R")
source("../lib/LoadKmertoneData.R")

# create reference table of kmers
sample.kmer <- data.table(kmer = GenerateKmerTable(k, ref = FALSE))
kmer.ref <- GenerateKmerTable(k, ref = TRUE)

# import kmertone results
import.kmertone <- function(action, file){
    out <- LoadKmertoneData(
        breakpoint.experiment = file,
        k = k,
        kmer.table = sample.kmer$kmer,
        kmer.ref.table = kmer.ref,
        action = action,
        to.target.dir = "../../02_Correlations/"
    )
    return(out)
}

actions <- c("z-score", "ratio")
out <- pbsapply(actions, function(action){
    out <- lapply(breakpoint.experiment, import.kmertone, action = action)
})

df.kmertone <- tibble(
    kmer = kmer.ref$kmer,
    z.score = out[[1]],
    ratio = out[[2]]
    ) %>% 
    rename_with(~c("kmer", "z.score", "ratio"))

N=10
df.subset <- rbind(
    df.kmertone %>% slice_max(order_by = z.score, n = N),
    df.kmertone %>% slice_min(order_by = z.score, n = N)
)

df.kmertone %>% 
    ggplot(aes(
        x = log2(ratio),
        y = abs(z.score)
    )) +
    geom_point(alpha = 0.4) + 
    ggrepel::geom_text_repel(
        data = df.subset,
        aes(
            x = log2(ratio),
            y = abs(z.score),
            label = kmer
        ),
        box.padding = 0.5,
        max.overlaps = Inf,
        min.segment.length = 0,
        size = 3
    ) + 
    labs(
        x = "Log2 probability ratios",
        y = "Absolute z-scores"
    )
