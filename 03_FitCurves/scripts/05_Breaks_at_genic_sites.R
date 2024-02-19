# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(cluster)))
suppressPackageStartupMessages(suppressWarnings(library(mclust)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults
pbapply::pboptions(char = "=")

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

#' 1. Import all genomic features for analysis
files_to_load <- list.files(
    path = "../../../DNAFragility/COSMIC/data/annotations/",
    pattern = "group_*",
    full.names = TRUE
)
group_names <- stringr::str_extract(
    string = files_to_load,
    pattern = "(?<=group_)[^.]+(?=\\.csv)"
)
all_groups <- lapply(files_to_load, fread, showProgress = FALSE)
names(all_groups) <- group_names

# only keep the below genic features
all_groups <- all_groups$genic_features
to_keep <- c(
    "CDS", "Exons", 
    "Five_prime_UTR", "Three_prime_UTR", 
    "Promoters", "TSS"
)
all_groups <- all_groups[(type %in% to_keep)]

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
        "../../data/00_Breakage/Ultrasonication/Simons_exp_1",
        "../../data/00_Breakage/Ancient_DNA/Altai_Neandertal",
        "../../data/00_Breakage/cfDNA/Ovarian_cancer",
        "../../data/00_Breakage/sBLISS/K562_Top2_mediated_DSBs/ETO",
        "../../data/00_Breakage/Enzymatic/EcoRV_HeLa_cells"
    ),
    hex = c(
        "#8B8000",
        "#9CDC9C",
        "#FFA435",
        "#BFABCB",
        "#BB357E"   
    )
)

#' 3. Calculate percent overlap with the genic region of interest
overlap_breaks <- lapply(1:nrow(df_rmsd), function(row){
    cat("Break class ", row, " of ", nrow(df_rmsd), "...\n")

    # get breakage class of interest
    t1 <- Sys.time()
    cur.msg <- paste0("Getting breakpoints and lifting over to hg38...")
    l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
    cat(paste0(cur.msg, l))

    df_bp <- lapply(1:22, function(chr){
        fetch.file <- paste0(
            df_rmsd$bp_exp[row], 
            "/breakpoint_positions/chr", 
            chr, ".csv"
        )
        df <- fread(file = fetch.file, showProgress = FALSE)
        if("freq" %in% colnames(df)) df[, lev.dist := NULL]
        setnames(df, "start")
        df[, `:=`(seqnames = paste0("chr", chr), width = 1)]
    })
    df_bp <- rbindlist(df_bp)
    df_bp <- plyranges::as_granges(df_bp)

    # genome annotations are in hg38, need to liftover breakpoints to hg38
    suppressPackageStartupMessages(suppressWarnings(library(rtracklayer)))
    chain <- import.chain("../data/liftover/hg19ToHg38.over.chain")
    df_bp <- liftOver(df_bp, chain)
    df_bp <- unlist(as(df_bp, "GRangesList"))

    total.time <- Sys.time() - t1
    cat("DONE! --", signif(total.time[[1]], 2), 
        attr(total.time, "units"), "\n")

    # get genic group
    types <- unique(all_groups$type)
    
    # calculate average overlap for all autosomes
    overlap_types <- lapply(1:length(types), function(x){
        t1 <- Sys.time()
        cur.msg <- paste0(
            "Calculating overlap with feature ", 
            x, " of ", length(types)
        )
        l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
        cat(paste0(cur.msg, l))

        df_temp <- plyranges::as_granges(all_groups[type == types[x]])

        overlaps <- GenomicRanges::findOverlaps(df_bp, df_temp)
        overlap_counts <- as.data.table(queryHits(overlaps))
        overlap_counts <- overlap_counts[, .(Total = .N), by = V1]
        setnames(overlap_counts, c("ID", "count"))
        overlaps <- sum(overlap_counts$count)/length(df_bp)*100

        overlaps <- data.table(
            "Genic_feature" = types[x], 
            "Perc_overlap" = overlaps
        )

        total.time <- Sys.time() - t1
        cat("DONE! --", signif(total.time[[1]], 2), 
            attr(total.time, "units"), "\n")

        gc()
        return(overlaps)
    })
    overlap_types <- rbindlist(overlap_types)
    overlap_types[, `:=`(
        Break_class = df_rmsd$class[row],
        hex = df_rmsd$hex[row]
    )]
    gc()
    return(overlap_types)
})
overlap_breaks <- rbindlist(overlap_breaks)

#' 4. Plot results 
p1 <- as_tibble(overlap_breaks) %>% 
    dplyr::mutate(
        Genic_feature = dplyr::case_when(
            Genic_feature == "Five_prime_UTR" ~ "5' UTR",
            Genic_feature == "Three_prime_UTR" ~ "3' UTR",
            TRUE ~ Genic_feature
        )
    ) %>% 
    dplyr::group_by(Genic_feature) %>% 
    dplyr::mutate(Median = median(Perc_overlap)) %>% 
    dplyr::arrange(desc(Median)) %>% 
    dplyr::mutate(Genic_feature = as.factor(Genic_feature)) %>% 
    dplyr::ungroup() %>% 
    ggplot(aes(x = forcats::fct_inorder(Genic_feature), y = Perc_overlap)) +
    geom_boxplot() + 
    geom_point(aes(color = hex), size = 3) + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    scale_color_identity() + 
    labs(
        x = "Genomic feature",
        y = "Overlap with breaks, %"
    )

ggsave(
    filename = "../figures/Breaks_Genomicfeatureoverlap.pdf",
    plot = p1,
    height = 6, 
    width = 6
)