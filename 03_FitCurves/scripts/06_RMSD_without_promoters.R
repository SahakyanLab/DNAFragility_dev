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
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults
pbapply::pboptions(char = "=")

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

source("../lib/SequenceEffect.R")
source("../lib/FitCurves.R")
Rcpp::sourceCpp("../lib/SeqEffect.cpp")

#' 1. Import all genomic features for analysis
files_to_load <- list.files(
    path = "../../../05_Cosmic/data/annotations",
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
    # "CDS", 
    # "Exons", 
    # "Five_prime_UTR", 
    # "Three_prime_UTR", 
    "Promoters"
    # "TSS"
)
all_groups <- all_groups[(type %in% to_keep)]
df_genic_feat <- plyranges::as_granges(all_groups)

# genome annotations are in hg38, need to liftover to match breakpoint's hg19
suppressPackageStartupMessages(suppressWarnings(library(rtracklayer)))
chain <- import.chain("../data/liftover/hg38ToHg19.over.chain")
df_genic_feat <- liftOver(df_genic_feat, chain)
df_genic_feat <- unlist(as(df_genic_feat, "GRangesList"))

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
    )
)

round_to_nearest_even <- function(x) round(x/2)*2
rmsd_range <- c(-501, 501)

#' 3. Calculate percent overlap with the genic region of interest
rmsd_results <- lapply(1:nrow(df_rmsd), function(row){
    t1 <- Sys.time()
    cat("Processing break class", row, "of", nrow(df_rmsd), "...\n")

    df_bp <- pbapply::pblapply(1:22, function(chr){
        fetch.file <- paste0(
            "../../data/",
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

    # remove direct overlap between break and promoter region
    df_bp_no_overlaps <- plyranges::filter_by_non_overlaps(
        df_bp, df_genic_feat
    ) %>% suppressWarnings()

    # remove overlaps witin the rmsd range too
    df_bp_no_overlaps_expand <- plyranges::stretch(
        plyranges::anchor_center(
            df_bp_no_overlaps %>% 
                dplyr::mutate(Breaks = start)
        ),
        abs(max(rmsd_range))
    )
    df_bp_no_overlaps_expand <- plyranges::filter_by_non_overlaps(
        df_bp_no_overlaps_expand, df_genic_feat
    ) %>% suppressWarnings()

    # get into the appropriate data format
    df_bp <- as.data.table(df_bp_no_overlaps_expand)
    df_bp[, `:=`(
        end = NULL, 
        width = NULL,
        strand = NULL,
        start = Breaks, 
        Breaks = NULL
    )]

    split_dfs <- split(df_bp, df_bp$seqnames)
    df_bp <- lapply(split_dfs, function(df) df[, seqnames := NULL])

    seq_effect <- SequenceEffect$new(
        chr = 1:22,
        return_vals = TRUE
    )
    seq_effect$calc_seq_effect(
        k = 8, 
        rmsd.range = rmsd_range,
        break_data = df_bp,
        assembly = "hg19"
    )

    # get average RMSD per position
    res <- lapply(1:length(seq_effect$results), function(chr){
        rmsd_vals <- seq_effect$results[[chr]]
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

    # fit Gaussian curves
    set.seed(1234)
    fit_curves <- FitCurves$new(
        chr = 1,
        rmsd_values = list("chr1" = res),
        return_vals = TRUE
    )
    fit_curves$generate_rmsd_plots(
        k = 8, 
        annot_plot = FALSE, 
        per_chr = TRUE
    )

    # map curves
    round_to_nearest_even <- function(x) round(x/2)*2

    map_ranges <- fit_curves$results$chr1$df %>% 
        dplyr::filter(rowid == "ranges") %>% 
        dplyr::select(-rowid)
    which_ind <- which(!is.na(map_ranges))
    map_ranges <- map_ranges[, which_ind]

    if(length(map_ranges) == 3){
        df_stats <- map_ranges %>% 
            tidyr::pivot_longer(
                cols = everything(),
                names_to = "Key",
                values_to = "Values"
            ) %>% 
            dplyr::mutate(
                Curves = dplyr::case_when(
                    Values == min(Values) ~ "short.range",
                    Values == median(Values) ~ "mid.range",
                    Values == max(Values) ~ "long.range"
                )
            ) %>% 
            dplyr::select(Key, Curves) %>% 
            tidyr::pivot_wider(
                names_from = Key,
                values_from = Curves
            )
    } else if(length(map_ranges) == 2){
        df_stats <- map_ranges %>% 
            tidyr::pivot_longer(
                cols = everything(),
                names_to = "Key",
                values_to = "Values"
            ) %>% 
            dplyr::mutate(
                Curves = dplyr::case_when(
                    Values == min(Values) ~ "short.range",
                    TRUE ~ NA_character_
                )
            )
            
        unused_range <- setdiff(
            c("short.range", "mid.range", "long.range"), 
            unique(df_stats$Curves)
        )

        # match with the closest one.
        map_to_closest_range <- fread(paste0(
            "../data/ranges/kmer_8_Ranges_",
            "cutoffs_from_clustering_all-exp.csv"
        ))
        map_to_closest_range <- map_to_closest_range[grepl(
            pattern = df_rmsd$bp_exp[row], 
            x = map_to_closest_range$exp
        ),]
        
        map_to_closest_range <- as_tibble(map_to_closest_range) %>% 
            dplyr::filter(rowid == "ranges") %>% 
            dplyr::select(-Curve, -exp, -rowid) %>% 
            tidyr::pivot_longer(
                cols = everything(),
                names_to = "Key",
                values_to = "Value"
            ) %>% 
            tidyr::drop_na() %>% 
            dplyr::mutate(
                Key = stringr::str_remove(
                    string = Key, pattern = "-1|-2"
                )
            ) %>% 
            dplyr::group_by(Key) %>% 
            dplyr::summarise(Value = max(Value)) %>% 
            dplyr::mutate(Value = round_to_nearest_even(Value)) %>% 
            dplyr::rename_with(~c('Curve', 'Cluster'))

        find_closest_curve <- function(value){
            return(
                map_to_closest_range %>% 
                    dplyr::mutate(Diff = abs(Cluster-value)) %>% 
                    dplyr::arrange(Diff) %>% 
                    dplyr::slice(1) %>% 
                    dplyr::pull(Curve)
            )
        }
        
        df_stats <- df_stats %>% 
            dplyr::mutate(
                Curves = ifelse(is.na(Curves), 
                    sapply(Values, find_closest_curve), 
                    Curves
                )
            ) %>%
            dplyr::select(Key, Curves) %>% 
            tidyr::pivot_wider(
                names_from = Key,
                values_from = Curves
            )
    } else {
        which_ind <- which(!is.na(map_ranges))[1]
        df_stats <- tibble("curve.one" = "short.range")
        colnames(df_stats) <- colnames(map_ranges)[which_ind]
    }

    df_stats <- setNames(
        as.character(df_stats),
        colnames(df_stats)
    )

    # get data tables into appropriate format
    output <- fit_curves$results$chr1
    limits <- nrow(output$all_gc)/2-1  
    x_vals <- (-limits):(nrow(output$all_gc)-limits-1)

    curves_vals <- data.table(output$curve_vals, x = x_vals)
    curves_vals[, col.label := as.character(df_stats[curves_vals$col.label])]

    output <- rbind(
        data.table(vals = res, col.label = "RMSD", x = x_vals),
        curves_vals,
        data.table(
            vals = output$all_gc$V1, 
            col.label = "full.range",
            x = x_vals
        )
    )
    setnames(output, c("rmsd", "curve", "x"))
    output[, `:=`(
        Break_class = df_rmsd$class[row],
        hex = df_rmsd$hex[row]
    )]

    # temp <- as_tibble(output) %>% 
    #     dplyr::mutate(
    #         hex = dplyr::case_when(
    #             curve == "RMSD" ~ "darkgrey",
    #             curve == "short.range" ~ "#619B61",
    #             curve == "mid.range" ~ "#355c7d",
    #             curve == "long.range" ~ "#BB357E",
    #             curve == "full.range" ~ "#f8b195"
    #         ),
    #         alpha = dplyr::case_when(
    #             curve == "RMSD" ~ 0.5,
    #             TRUE ~ 1
    #         ),
    #         lw = dplyr::case_when(
    #             curve == "RMSD" ~ 0.5,
    #             curve == "short.range" ~ 1,
    #             curve == "mid.range" ~ 1,
    #             curve == "long.range" ~ 1,
    #             curve == "full.range" ~ 1.2
    #         ),
    #         curve = factor(curve, levels = c(
    #             "RMSD",
    #             "curve.one", 
    #             "mid.range", 
    #             "long.range", 
    #             "full.range"
    #             )
    #         )
    #     ) 

    # p1 <- temp %>% 
    #     ggplot(aes(x = x, y = rmsd, col = hex, alpha = alpha)) +
    #     geom_line(linewidth = temp$lw) +
    #     scale_color_identity() +
    #     theme_bw() + 
    #     theme_classic() 
    # ggsave(
    #     filename = "../figures/demo.pdf",
    #     plot = p1,
    #     height = 6, 
    #     width = 6
    # )
    return(output)
})
rmsd_results <- rbindlist(rmsd_results)

dat <- as_tibble(rmsd_results) %>% 
    dplyr::mutate(
        hex = dplyr::case_when(
            curve == "RMSD" ~ "darkgrey",
            curve == "short.range" ~ "#619B61",
            curve == "mid.range" ~ "#355c7d",
            curve == "long.range" ~ "#BB357E",
            curve == "full.range" ~ "#f8b195"
        ),
        alpha = dplyr::case_when(
            curve == "RMSD" ~ 0.5,
            TRUE ~ 1
        ),
        lw = dplyr::case_when(
            curve == "RMSD" ~ 0.5,
            curve == "short.range" ~ 1,
            curve == "mid.range" ~ 1,
            curve == "long.range" ~ 1,
            curve == "full.range" ~ 1.2
        ),
        curve = factor(curve, levels = c(
            "RMSD",
            "curve.one", 
            "mid.range", 
            "long.range", 
            "full.range"
            )
        )
    )

# plot results
p1 <- dat %>%
    dplyr::mutate(
        Break_class = forcats::fct_inorder(Break_class)
    ) %>% 
    ggplot(aes(x = x, y = rmsd, col = hex)) + 
    geom_line(aes(alpha = alpha), linewidth = dat$lw) +
    facet_wrap(
        vars(Break_class), 
        ncol = 2,
        scales = "free"
    ) +
    scale_color_identity() + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        legend.position = "none"
    ) + 
    labs(
        x = "Position away from breakpoint, bp",
        y = "RMSD"
    )

ggsave(
    filename = "../figures/RMSD_without_promoters.pdf",
    plot = p1,
    height = 12, 
    width = 8
)