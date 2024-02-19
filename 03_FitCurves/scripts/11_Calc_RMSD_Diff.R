# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults
pbapply::pboptions(char = "=")

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

#' 2. Import RMSD values
df_rmsd <- tibble(
    class = c(
        "Mechanical",
        "Natural decay",
        "Cell free DNA",
        "Endogenous DSBs",
        "Enzymatic"
    ),
    bp_exp = c(
        "00_Breakage/Ultrasonication/Simons_exp_1",
        "00_Breakage/Ancient_DNA/Altai_Neandertal",
        "00_Breakage/cfDNA/Ovarian_cancer",
        "00_Breakage/sBLISS/K562_Top2_mediated_DSBs/ETO",
        "00_Breakage/Enzymatic/EcoRV_HeLa_cells"
    ),
    hex = c(
        "#8B8000",
        "#9CDC9C",
        "#FFA435",
        "#BFABCB",
        "#BB357E"   
    ),
    ref = c(
        "hs37d5", 
        "hg19",
        "hs37d5",
        "hg19",
        "hg19"
    )
)

row <- 1

kmer <- 8
rmsd_range <- -501:501
all_chr <- 1:22
df_rmsd <- df_rmsd[row,]
round_to_nearest_even <- function(x) round(x/2)*2
half_width <- max(abs(max(rmsd_range)))+ceiling(kmer/2)

# get the start and end positions of each chromosome
temp <- switch(df_rmsd$ref,
    "hg18" = BSgenome.Hsapiens.UCSC.hg18::BSgenome.Hsapiens.UCSC.hg18,
    "hg19" = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    "hg38" = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    "hs37d5" = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
)
refseq.table <- as.data.frame(temp@seqinfo)
refseq.table <- refseq.table[grepl(
    pattern = "^chr([1-9]|1[0-9]|2[0-2])$", 
    x = rownames(refseq.table)
),]

if(nrow(refseq.table) == 0){
    refseq.table <- as.data.frame(temp@seqinfo)
    refseq.table <- refseq.table[grepl(
        pattern = "^([1-9]|1[0-9]|2[0-2])$", 
        x = rownames(refseq.table)
    ),]
    rownames(refseq.table) <- paste0("chr", rownames(refseq.table))
}
chrs_len <- refseq.table$seqlengths         
chrs_len <- setNames(chrs_len, rownames(refseq.table))

# import breakpoint data
df_bp <- pbapply::pblapply(all_chr, function(chr){
    fetch.file <- paste0(
        "../../data/",
        df_rmsd$bp_exp[row], 
        "/breakpoint_positions/chr", 
        chr, ".csv"
    )
    df <- fread(file = fetch.file, showProgress = FALSE)
    if("freq" %in% colnames(df)) df[, lev.dist := NULL]
    setnames(df, "start")
    
    # make sure no position will be out-of-bounds
    df[, `:=`(
        lower_end = start-half_width,
        upper_end = start+half_width
    )]
    df[, to_keep := ifelse(
        lower_end > 0 & upper_end < chrs_len[chr],
        TRUE, FALSE
    )]
    df <- df[to_keep == TRUE]
    df[, `:=`(lower_end = NULL, upper_end = NULL, to_keep = NULL)]
})
names(df_bp) <- paste0("chr", all_chr)
    
RMSD_method_all <- c(
    "right_vs_left",
    # "right_vs_mid",
    "right_and_left"
)

path_to_save <- "../data/ranges/demo"
dir.create(path = path_to_save, showWarnings = FALSE)

results <- lapply(RMSD_method_all, function(RMSD_method){
    if(RMSD_method == "right_vs_left"){
        Rcpp::sourceCpp("../lib/SeqEffect.cpp")
    } else if(RMSD_method == "right_vs_mid"){
        Rcpp::sourceCpp("../lib/SeqEffect_DiffToMid.cpp")
    } else if(RMSD_method == "right_and_left"){
        Rcpp::sourceCpp("../lib/SeqEffect_DiffRightAndLeft.cpp")
    }

    # to save RMSD results
    results <- pbapply::pblapply(all_chr, function(chr){
        ref_rmsd <- paste0(
            "../../data/ref/", df_rmsd$ref[row], 
            "/chr", chr, ".fasta.gz"
        )

        return(calc_kmer_freq(
            bp_pos = df_bp[[paste0("chr", chr)]]$start,
            filename = ref_rmsd,
            kmer = kmer,
            rmsd_range = rmsd_range
        ))
    })
    names(results) <- paste0("chr", all_chr)

    # get average RMSD per position
    res <- lapply(all_chr, function(chr){
        rmsd_vals <- results[[chr]]
        limits <- length(rmsd_vals)/2-1  
        temp <- data.table(
            seqnames = paste0("chr", chr),
            x = (-limits):(length(rmsd_vals)-limits-1),
            vals = rmsd_vals
        )
        return(temp)
    })
    res <- rbindlist(res)
    res <- as_tibble(res) %>% 
        dplyr::group_by(x) %>% 
        dplyr::summarise(y = mean(vals, na.rm = TRUE)) %>% 
        dplyr::ungroup() %>%
        dplyr::pull(y)

    # get data tables into appropriate format
    limits <- length(res)/2-1  
    x_vals <- (-limits):(length(res)-limits-1)
    output <- data.table(
        x = x_vals,
        rmsd = res,
        Break_class = df_rmsd$class[row],
        hex = df_rmsd$hex[row],
        RMSD_method = RMSD_method
    )
    fwrite(output, paste0(path_to_save, "/RMSD_Demo_", RMSD_method, ".csv"))

    p1 <- as_tibble(output) %>% 
        ggplot(aes(x = x, y = rmsd)) + 
        geom_line(linewidth = 0.8) + 
        coord_cartesian(ylim = c(0, NA)) + 
        theme_bw() + 
        theme_classic() + 
        theme(text = element_text(size = 15)) + 
        labs(
            x = "Position away from breakpoint, bp",
            y = "RMSD",
            subtitle = paste0("RMSD method: ", RMSD_method)
        )

    ggsave(
        filename = paste0(path_to_save, "/RMSD_Demo_", RMSD_method, ".pdf"),
        plot = p1,
        height = 7, width = 8
    )

    return(output)
})
results <- rbindlist(results)

p1 <- as_tibble(results) %>% 
    dplyr::mutate(
        RMSD_method = ifelse(
            RMSD_method == "right_and_left", 
            "i - [i-1] and i - [i+1]",
            "i - [i-1] or i - [i+1]"
        )
    ) %>% 
    ggplot(aes(x = x, y = rmsd)) + 
    geom_line(linewidth = 0.8) + 
    # coord_cartesian(ylim = c(0, NA)) + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    facet_wrap(vars(RMSD_method), nrow = 1) + 
    scale_y_continuous(labels = scales::label_number(scale = 1e6)) +
    labs(
        x = "Position away from breakpoint, bp",
        y = expression("RMSD, x10"^-6*"")
    )

ggsave(
    filename = paste0(path_to_save, "/RMSD_Demo_All.pdf"),
    plot = p1,
    height = 5, width = 10
)