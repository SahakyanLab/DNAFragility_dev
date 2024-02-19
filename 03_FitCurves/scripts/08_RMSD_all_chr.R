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

#' 1. Import RMSD values
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
    ind = c(
        1,
        28,
        81,
        56,
        35
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
chr <- 1:22

#' 3. Calculate percent overlap with the genic region of interest
rmsd_results <- lapply(1:nrow(df_rmsd), function(row){
    t1 <- Sys.time()
    cat(paste0("Processing break class ", row, " of ", nrow(df_rmsd), "...\n"))

    seq_effect <- SequenceEffect$new(
        chr = chr,
        which_exp_ind = df_rmsd$ind[row],
        return_vals = TRUE
    )
    seq_effect$calc_seq_effect(
        k = 8, 
        rmsd.range = c(-501, 501)
    )   

    # fit Gaussian curves
    set.seed(1234)
    fit_curves <- FitCurves$new(
        chr = chr,
        rmsd_values = seq_effect$results,
        return_vals = TRUE
    )
    fit_curves$generate_rmsd_plots(k = 8, per_chr = TRUE)

    results <- lapply(1:length(fit_curves$results), function(x){
        # map curves
        round_to_nearest_even <- function(x) round(x/2)*2

        map_ranges <- fit_curves$results[[x]]$df %>% 
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
                pattern = stringr::str_remove(
                    string = df_rmsd$bp_exp[row],
                    pattern = "../../data/"
                ),
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
        
        # rmsd curve values
        rmsd_vals <- seq_effect$results[[x]]
        limits <- length(rmsd_vals)/2-1

        dat <- fit_curves$results[[x]]$curve_vals
        nr.of.curves <- length(unique(dat$col.label))
        x_vals <- -limits:(length(rmsd_vals)-limits-1)
        x_vals_rep <- rep(x_vals, nr.of.curves)
        dat[, x := x_vals_rep]
        setnames(dat, c("y", "curve", "x"))

        dat[, curve := as.character(df_stats[dat$curve])]
        dat <- rbind(
            data.table(
                y = rmsd_vals, 
                curve = "RMSD", 
                x = x_vals
            ),
            dat, 
            data.table(
                y = fit_curves$results[[x]]$all_gc$V1, 
                curve = "full.range", 
                x = x_vals
            )
        )
        dat[, Break_class := df_rmsd$class[row]]
        dat$chr <- chr[x]

        return(dat)
    })
    df_vals <- rbindlist(results)

    # smalltest val
    min_num <- min(df_vals$y, na.rm = TRUE)
    exponent <- floor(log10(abs(min_num)))+2
    scaling_fac <- 10^-(exponent)
    scale_labels <- function(x) sprintf("%.1f", x * scaling_fac)

    # plot results
    plot_rmsd <- as_tibble(df_vals) %>% 
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
                curve == "full.range" ~ 1.5
            ),
            curve = factor(curve, levels = c(
                    "RMSD",
                    "short.range", 
                    "mid.range", 
                    "long.range", 
                    "full.range"
                )
            ),
            chr = paste0("chr", chr),
            chr = factor(chr, levels = c(paste0("chr", 1:22)))
        ) %>% 
        ggplot(aes(x = x, y = y, col = hex, alpha = alpha)) + 
        geom_line(linewidth = 1.2) +
        facet_wrap(vars(chr), ncol = 6, scales = "free_y") +
        scale_color_identity() + 
        scale_y_continuous(labels = scale_labels) +
        theme_bw() + 
        theme_classic() +
        theme(
            text = element_text(size = 15),
            legend.position = "none"
        ) + 
        labs(
            x = "Position away from breakpoint, bp",
            y = bquote("RMSD, x10"^ .(exponent))
        )

    ggsave(
        filename = paste0(
            "../figures/RMSD_all_chr_", 
            df_rmsd$class[row], ".pdf"
        ),
        plot = plot_rmsd,
        height = 8, 
        width = 16
    )

    return(df_vals)
})
df_vals <- rbindlist(rmsd_results)

# save values
fwrite(
    df_vals,
    file = "../data/RMSD_all_chr_selected_example_datasets.csv"
)