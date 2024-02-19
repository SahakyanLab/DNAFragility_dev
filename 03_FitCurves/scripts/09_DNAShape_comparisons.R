# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(DNAshapeR)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults
pbapply::pboptions(char = "=")

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

source("../lib/SequenceEffect.R")
source("../lib/FitCurves.R")
Rcpp::sourceCpp("../lib/SeqEffect.cpp")

# obtain A/B subcompartments
file_names <- c(
    "../../data/HiC_annotations/K562_track_hg19.bed",
    "../../data/HiC_annotations/HeLa_track_hg19.bed"
)
csv_files <- gsub("bed", "csv", file_names)
if(!all(file.exists(csv_files))){
    for(f in file_names){
        df_comp <- plyranges::read_bed(f)
        df_comp <- as_tibble(df_comp) %>%
            dplyr::select(seqnames, start, end, name) %>% 
            dplyr::mutate(
                name = dplyr::case_when(
                    grepl("^A[0-9]", name) ~ "open",
                    grepl("^B[0-9]", name) ~ "closed"
                )
            ) %>% 
            dplyr::group_split(name)

        for(comp in 1:length(df_comp)){
            temp <- df_comp[[comp]]
            temp <- as.data.table(temp)
            which_comp <- unique(temp$name)
            temp[, name := NULL]
            fwrite(
                temp,
                gsub(".bed", paste0("_", which_comp, ".csv"), f)
            )
        }
    }
}

# # main study
# df <- tibble(
#     class = c(
#         "EcoRV_full", "EcoRV_A_compartment", "EcoRV_B_compartment", 
#         "Nt_BbvCI_full", "Nt_BbvCI_A_compartment", "Nt_BbvCI_B_compartment"
#     ),
#     break_type = rep("enzymatic", 6),
#     bp_exp = c(
#         rep("../../data/00_Breakage/Enzymatic/EcoRV_HeLa_cells", 3),
#         rep("../../data/00_Breakage/Enzymatic/Nt_BbvCI_K562_cells/Sap", 3)
#     ),
#     compartments = c(
#         "",
#         "../../data/HiC_annotations/HeLa_track_hg19_open.csv",
#         "../../data/HiC_annotations/HeLa_track_hg19_closed.csv",
#         "",
#         "../../data/HiC_annotations/K562_track_hg19_open.csv",
#         "../../data/HiC_annotations/K562_track_hg19_closed.csv"
#     ),
#     ranges = c(
#         rep("../data/00_Breakage/Enzymatic/EcoRV_HeLa_cells/key_stats_kmer_8.csv", 3),
#         rep("../data/00_Breakage/Enzymatic/Nt_BbvCI_K562_cells/Sap/key_stats_kmer_8.csv", 3)
#     ),
#     ind = c(rep(35, 3), rep(37, 3)),
#     hex = rep("#BB357E", 6)
# )

# main study
df <- tibble(
    class = c("EcoRV", "Nt_BbvCI"),
    break_type = "enzymatic",
    bp_exp = c(
        "../../data/00_Breakage/Enzymatic/EcoRV_HeLa_cells",
        "../../data/00_Breakage/Enzymatic/Nt_BbvCI_K562_cells/Sap"
    ),
    compartments = c("", ""),
    ranges = c(
        "../data/00_Breakage/Enzymatic/EcoRV_HeLa_cells/key_stats_kmer_8.csv", 
        "../data/00_Breakage/Enzymatic/Nt_BbvCI_K562_cells/Sap/key_stats_kmer_8.csv"
    ),
    ind = c(35, 37),
    hex = "#BB357E"
)

round_to_nearest_even <- function(x) round(x/2)*2
hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
seed <- 1234
mid_range <- 300
rmsd_range <- c(-501, 501)

# Run DNA shape calculation
rmsd_stats <- lapply(1:nrow(df), function(row){
    t1 <- Sys.time()
    cat(paste0("Processing break class ", row, " of ", nrow(df), "...\n"))

    # import breakpoint data
    df_bp <- lapply(1:22, function(chr){
        df_bp <- fread(paste0(
            df$bp_exp[row], 
            "/breakpoint_positions/chr", 
            chr, ".csv"
        ))
        setnames(df_bp, "start")
        df_bp[, `:=`(
            seqnames = paste0("chr", chr),
            width = 1,
            strand = "+"
        )]
        return(df_bp)
    })
    df_bp <- plyranges::as_granges(rbindlist(df_bp))

    # filter for compartment of interest 
    if(df$compartments[row] != ""){
        df_compartment <- fread(df$compartments[row])
        df_compartment <- plyranges::as_granges(df_compartment)
        df_bp <- plyranges::filter_by_overlaps(df_bp, df_compartment)
    }

    # extract sequences
    df_bp_MR <- plyranges::stretch(
        plyranges::anchor_center(df_bp),
        mid_range-1
    )

    # pass through DNAShapeR calculation
    DNAshapeR::getFasta(
        df_bp_MR, hg19, 
        width = as.numeric(width(df_bp_MR)[1]),
        filename = "temp.fa"
    ) %>% suppressMessages()
    pred <- DNAshapeR::getShape("temp.fa") %>% suppressMessages()

    # compute summary statistics
    shape_name <- c("HelT", "MGW", "ProT", "Roll")
    shape_plots <- lapply(1:length(shape_name), function(i){
        dat <- pred[[shape_name[i]]]

        # remove first and second columns as they are NAs
        dat <- dat[, -c(1, ncol(dat))]
        dat <- as.data.frame(dat)
        limits <- ceiling(ncol(dat)/2-1)
        colnames(dat) <- as.character(-limits:(ncol(dat)-limits-1))

        # Calculate the average and confidence interval for each column
        dat_summary <- dat %>% 
            dplyr::summarise(across(
                everything(),
                list(
                    Mean = ~ mean(.x, na.rm = TRUE),
                    Lower_CI = ~ quantile(.x, probs = 0.025, na.rm = TRUE),
                    Upper_CI = ~ quantile(.x, probs = 0.975, na.rm = TRUE)
                )
            ))

        dat_long <- dat_summary %>% 
            tidyr::pivot_longer(
                cols = everything(), 
                names_to = c("x_val", "Column"), 
                names_pattern = "(-?\\d+)_(Mean|Lower_CI|Upper_CI)"
            ) %>% 
            dplyr::mutate(x_val = as.numeric(x_val))

        dat_wide <- dat_long %>% 
            tidyr::pivot_wider(
                names_from = "Column",
                values_from = "value"
            ) %>% 
            dplyr::mutate(shape = shape_name[i])

        return(dat_wide)
    })
    shape_plots <- do.call(dplyr::bind_rows, shape_plots)

    dir.create(
        "../figures/DNA_shape_experiments",
        showWarnings = FALSE,
        recursive = TRUE
    )
    plot_shapes <- shape_plots %>% 
        ggplot(aes(x = x_val)) +
        geom_ribbon(
            aes(ymin = Lower_CI, ymax = Upper_CI),
            fill = "grey80",
            alpha = 0.5
        ) + 
        geom_line(
            aes(y = Mean),
            color = "#B33F40",
            linewidth = 1.2
        ) + 
        facet_wrap(vars(shape), scales = "free_y", ncol = 2) + 
        theme_bw() + 
        theme_classic() + 
        theme(text = element_text(size = 20)) + 
        labs(
            x = "Position away from breakpoint, bp",
            y = "Value of shape"
        )

    ggsave(
        filename = paste0(
            "../figures/DNA_shape_experiments/",
            df$class[row], "_shapes.pdf"
        ),
        plot = plot_shapes,
        height = 10, width = 14
    )

    # GC content
    df_bp_MR_seq <- getSeq(hg19, df_bp_MR)
    GC_mat <- pbapply::pbsapply(1:length(df_bp_MR_seq), function(x){
        letter.mat <- letterFrequencyInSlidingView(
            x = df_bp_MR_seq[[x]],
            view.width = 1, 
            letters="GC", 
            OR=0
        )
        return(matrixStats::rowSums2(letter.mat, na.rm = TRUE))
    })
    GC_mat <- t(GC_mat)
    GC_mat_colsum <- matrixStats::colSums2(GC_mat)/length(df_bp_MR_seq)
    limits <- ceiling(length(GC_mat_colsum)/2-1)
    df_GC <- tibble(
            x = -limits:(length(GC_mat_colsum)-limits-1),
            GC = GC_mat_colsum * 100
        ) 
        
    plot_GC <- df_GC %>% 
        ggplot(aes(x = x, y = GC)) + 
        geom_line(color = "#B33F40", linewidth = 1.2) + 
        theme_bw() + 
        theme_classic() + 
        theme(text = element_text(size = 15)) + 
        labs(
            x = "Position away from breakpoint, bp",
            y = "GC content, %"
        )

    ggsave(
        filename = paste0(
            "../figures/DNA_shape_experiments/",
            df$class[row], "_GC.pdf"
        ),
        plot = plot_GC,
        height = 6, width = 8
    )

    # # RMSD curves
    # df_bp <- as.data.table(df_bp)
    # df_bp[, `:=`(end = NULL, width = NULL, strand = NULL)]

    # split_dfs <- split(df_bp, df_bp$seqnames)
    # df_bp <- lapply(split_dfs, function(df) df[, seqnames := NULL])

    # seq_effect <- SequenceEffect$new(
    #     chr = 1:22,
    #     return_vals = TRUE
    # )
    # seq_effect$calc_seq_effect(
    #     k = 8, 
    #     rmsd.range = rmsd_range,
    #     break_data = df_bp,
    #     assembly = "hg19"
    # )

    # # get average RMSD per position
    # res <- lapply(1:length(seq_effect$results), function(chr){
    #     rmsd_vals <- seq_effect$results[[chr]]
    #     limits <- length(rmsd_vals)/2-1  
    #     temp <- data.table(
    #         seqnames = paste0("chr", chr),
    #         x = (-limits):(length(rmsd_vals)-limits-1),
    #         vals = rmsd_vals
    #     )
    #     return(temp)
    # })
    # res <- rbindlist(res)
    # res <- as_tibble(res) %>% 
    #     dplyr::group_by(x) %>% 
    #     dplyr::summarise(y = mean(vals, na.rm = TRUE)) %>% 
    #     dplyr::ungroup() %>%
    #     dplyr::pull(y)

    # # fit Gaussian curves
    # set.seed(seed)
    # fit_curves <- FitCurves$new(
    #     chr = 1,
    #     rmsd_values = list("chr1" = res),
    #     return_vals = TRUE
    # )
    # fit_curves$generate_rmsd_plots(
    #     k = 8, 
    #     annot_plot = FALSE, 
    #     per_chr = TRUE
    # )

    # output <- fit_curves$results$chr1
    # limits <- nrow(output$all_gc)/2-1  
    # x_vals <- (-limits):(nrow(output$all_gc)-limits-1)
    # output <- rbind(
    #     data.table(vals = res, col.label = "RMSD", x = x_vals),
    #     data.table(output$curve_vals, x = x_vals),
    #     data.table(
    #         vals = output$all_gc$V1, 
    #         col.label = "curve.all",
    #         x = x_vals
    #     )
    # )
    # setnames(output, c("rmsd", "curve", "x"))
    # output[, `:=`(
    #     Break_class = df$class[row],
    #     hex = df$hex[row]
    # )]

    # dat <- as_tibble(output) %>% 
    #     dplyr::mutate(
    #         hex = dplyr::case_when(
    #             curve == "RMSD" ~ "darkgrey",
    #             curve == "curve.one" ~ "#619B61",
    #             curve == "curve.two" ~ "#355c7d",
    #             curve == "curve.three" ~ "#BB357E",
    #             curve == "curve.all" ~ "#f8b195"
    #         ),
    #         alpha = dplyr::case_when(
    #             curve == "RMSD" ~ 0.5,
    #             TRUE ~ 1
    #         ),
    #         lw = dplyr::case_when(
    #             curve == "RMSD" ~ 0.5,
    #             curve == "curve.one" ~ 1,
    #             curve == "curve.two" ~ 1,
    #             curve == "curve.three" ~ 1,
    #             curve == "curve.all" ~ 1.2
    #         ),
    #         curve = factor(curve, levels = c(
    #             "RMSD",
    #             "curve.one", 
    #             "curve.two", 
    #             "curve.three", 
    #             "curve.all"
    #             )
    #         )
    #     )
        
    # # plot results
    # p1 <- dat %>%
    #     dplyr::mutate(
    #         Break_class = forcats::fct_inorder(Break_class)
    #     ) %>% 
    #     ggplot(aes(x = x, y = rmsd, col = hex)) + 
    #     geom_line(aes(alpha = alpha), linewidth = dat$lw) +
    #     facet_wrap(
    #         vars(Break_class), 
    #         ncol = 2,
    #         scales = "free"
    #     ) +
    #     scale_color_identity() + 
    #     theme_bw() + 
    #     theme_classic() + 
    #     theme(
    #         text = element_text(size = 15),
    #         legend.position = "none"
    #     ) + 
    #     labs(
    #         x = "Position away from breakpoint, bp",
    #         y = "RMSD"
    #     )

    # ggsave(
    #     filename = paste0(
    #         "../figures/DNA_shape_experiments/",
    #         df$class[row], "_RMSD.pdf"
    #     ),
    #     plot = p1,
    #     height = 5, 
    #     width = 6
    # )
    
    # # curve stats
    # curve_stats <- as.data.table(fit_curves$results$chr1$df)
    # setnames(curve_stats, c("ID", "Curve_one", "Curve_two", "Curve_three"))
    # curve_stats[, Break_class := df$class[row]]

    # df_GC <- df_GC %>% 
    #     dplyr::mutate(
    #         Break_class = df$class[row],
    #         hex = df$hex[row]
    #     )

    # shape_plots <- shape_plots %>% 
    #     dplyr::mutate(
    #         Break_class = df$class[row],
    #         hex = df$hex[row]
    #     )

    # return(list(curve_stats, output, df_GC, shape_plots))
    # return(list(curve_stats, output))
})

curve_stats <- rbindlist(sapply(rmsd_stats, `[`, 1))
df_rmsd <- rbindlist(sapply(rmsd_stats, `[`, 2))
# df_GC <- do.call(rbind, sapply(rmsd_stats, `[`, 3))
# df_shape <- do.call(rbind, sapply(rmsd_stats, `[`, 4))

# RMSD plots
df_rmsd <- as_tibble(df_rmsd) %>% 
    dplyr::mutate(
        hex = dplyr::case_when(
            curve == "RMSD" ~ "darkgrey",
            curve == "curve.one" ~ "#619B61",
            curve == "curve.two" ~ "#355c7d",
            curve == "curve.three" ~ "#BB357E",
            curve == "curve.all" ~ "#f8b195"
        ),
        alpha = dplyr::case_when(
            curve == "RMSD" ~ 0.5,
            TRUE ~ 1
        ),
        lw = dplyr::case_when(
            curve == "RMSD" ~ 0.5,
            curve == "curve.one" ~ 1,
            curve == "curve.two" ~ 1,
            curve == "curve.three" ~ 1,
            curve == "curve.all" ~ 1.2
        ),
        curve = factor(curve, levels = c(
            "RMSD",
            "curve.one", 
            "curve.two", 
            "curve.three", 
            "curve.all"
            )
        )
    )
    
# plot results
plot_rmsd <- df_rmsd %>%
    dplyr::mutate(
        Break_class = factor(Break_class, levels = c(
            "EcoRV_full", "Nt_BbvCI_full", 
            "EcoRV_A_compartment", "Nt_BbvCI_A_compartment",
            "EcoRV_B_compartment", "Nt_BbvCI_B_compartment"
        ))
    ) %>% 
    ggplot(aes(x = x, y = rmsd, col = hex)) + 
    geom_line(aes(alpha = alpha), linewidth = 1.3) +
    facet_wrap(
        vars(Break_class), 
        ncol = 2,
        scales = "free_y"
    ) +
    scale_color_identity() + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 20),
        legend.position = "none"
    ) + 
    labs(
        x = "Position away from breakpoint, bp",
        y = "RMSD"
    )

ggsave(
    filename = paste0(
        "../figures/DNA_shape_experiments/All_RMSD.pdf"
    ),
    plot = plot_rmsd,
    height = 12, width = 10
)

# # GC content plots
# plot_GC <- df_GC %>% 
#     dplyr::mutate(
#         Break_class = factor(Break_class, levels = c(
#             "EcoRV_full", "Nt_BbvCI_full", 
#             "EcoRV_A_compartment", "Nt_BbvCI_A_compartment",
#             "EcoRV_B_compartment", "Nt_BbvCI_B_compartment"
#         ))
#     ) %>% 
#     ggplot(aes(x = x, y = GC)) + 
#     geom_line(color = "#B33F40", linewidth = 1.2) + 
#     facet_wrap(vars(Break_class), ncol = 2) + 
#     theme_bw() + 
#     theme_classic() + 
#     theme(text = element_text(size = 20)) + 
#     coord_cartesian(ylim = c(0, 100)) + 
#     labs(
#         x = "Position away from breakpoint, bp",
#         y = "GC content, %"
#     )

# ggsave(
#     filename = paste0(
#         "../figures/DNA_shape_experiments/All_GC.pdf"
#     ),
#     plot = plot_GC,
#     height = 12, width = 10
# )

# GC content plots
plot_GC <- df_GC %>% 
    dplyr::filter(grepl("EcoRV_full|Nt_BbvCI_full", Break_class)) %>% 
    dplyr::mutate(
        Break_class = factor(Break_class, levels = c(
            "EcoRV_full", "Nt_BbvCI_full"
        ))
    ) %>% 
    ggplot(aes(x = x, y = GC)) + 
    geom_line(color = "#B33F40", linewidth = 1.2) + 
    facet_wrap(vars(Break_class), ncol = 2) + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 20)) + 
    coord_cartesian(ylim = c(0, 100)) + 
    labs(
        x = "Position away from breakpoint, bp",
        y = "GC content, %"
    )

ggsave(
    filename = paste0(
        "../figures/DNA_shape_experiments/Full_GC.pdf"
    ),
    plot = plot_GC,
    height = 6, width = 6
)

# # DNA shape plots
# plot_shapes <- df_shape %>% 
#     dplyr::mutate(
#         Break_class = factor(Break_class, levels = c(
#             "EcoRV_full", "Nt_BbvCI_full", 
#             "EcoRV_A_compartment", "Nt_BbvCI_A_compartment",
#             "EcoRV_B_compartment", "Nt_BbvCI_B_compartment"
#         ))
#     ) %>% 
#     ggplot(aes(x = x_val)) +
#     geom_ribbon(
#         aes(ymin = Lower_CI, ymax = Upper_CI),
#         fill = "grey80",
#         alpha = 0.5
#     ) + 
#     geom_line(
#         aes(y = Mean),
#         color = "#B33F40",
#         linewidth = 1.2
#     ) + 
#     facet_grid(
#         cols = vars(Break_class),
#         rows = vars(shape), 
#         scales = "free_y"
#     ) + 
#     theme_bw() + 
#     theme_classic() + 
#     theme(text = element_text(size = 20)) + 
#     labs(
#         x = "Position away from breakpoint, bp",
#         y = ""
#     )

# ggsave(
#     filename = paste0(
#         "../figures/DNA_shape_experiments/All_shapes.pdf"
#     ),
#     plot = plot_shapes,
#     height = 10, width = 19
# )

# DNA shape plots
plot_shapes <- df_shape %>% 
    dplyr::filter(grepl("EcoRV_full|Nt_BbvCI_full", Break_class)) %>% 
    dplyr::mutate(
        Break_class = factor(Break_class, levels = c(
            "EcoRV_full", "Nt_BbvCI_full"
        ))
    ) %>% 
    ggplot(aes(x = x_val)) +
    geom_ribbon(
        aes(ymin = Lower_CI, ymax = Upper_CI),
        fill = "grey80",
        alpha = 0.5
    ) + 
    geom_line(
        aes(y = Mean),
        color = "#B33F40",
        linewidth = 1.2
    ) + 
    facet_grid(
        cols = vars(Break_class),
        rows = vars(shape), 
        scales = "free_y"
    ) + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 20)) + 
    labs(
        x = "Position away from breakpoint, bp",
        y = ""
    )

ggsave(
    filename = paste0(
        "../figures/DNA_shape_experiments/Full_shapes.pdf"
    ),
    plot = plot_shapes,
    height = 10, width = 6
)