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

source("../lib/SequenceEffect.R")
source("../lib/FitCurves.R")
Rcpp::sourceCpp("../lib/SeqEffect.cpp")

set.seed(1234)
rng <- sample(100000, size = 10, replace = FALSE)
chr <- 1:22

results <- lapply(1:length(rng), function(x){
    seq_effect <- SequenceEffect$new(
        chr = chr,
        which_exp_ind = NULL,
        control = TRUE,
        return_vals = TRUE,
        seed = rng[x]
    )
    seq_effect$calc_seq_effect(
        k = 8, 
        rmsd.range = c(-501, 501),
        from_file = FALSE
    )
    df_bp <- seq_effect$results

    # get average RMSD per position
    res <- lapply(names(df_bp), function(chr){
        rmsd_vals <- df_bp[[chr]]
        limits <- length(rmsd_vals)/2-1  
        temp <- data.table(
            seqnames = chr,
            x = (-limits):(length(rmsd_vals)-limits-1),
            vals = rmsd_vals,
            size = seq_effect$size_of_data[[chr]]
        )
        return(temp)
    })
    res <- rbindlist(res)
    setnames(res, c("seqnames", "x", "vals", "size"))
    res <- as_tibble(res) %>% 
        dplyr::group_by(x) %>% 
        dplyr::summarise(
            y = mean(vals, na.rm = TRUE),
            size = ceiling(mean(size))
        ) %>% 
        dplyr::ungroup()

    limits <- nrow(res)/2-1
    results <- data.table(
        x = -limits:(nrow(res)-limits-1),
        rmsd_values = res$y,
        study_ID = paste0("ID_", rng[x]),
        size = res$size
    )

    p1 <- as_tibble(results) %>% 
        dplyr::mutate(rmsd_values = rmsd_values / sum(size)) %>%
        ggplot(aes(x = x, y = rmsd_values, col = 'darkred')) + 
        geom_line(linewidth = 0.8) + 
        scale_color_identity() + 
        theme_bw() + 
        theme_classic() + 
        labs(
            x = "Position away from breakpoint, bp",
            y = "RMSD"
        )
    ggsave(filename = "../figures/quick_demo.pdf", plot = p1)

    return(results)
})
results <- rbindlist(results)

#' 2. Plot results
p1 <- as_tibble(results) %>% 
    dplyr::mutate(
        rmsd_values = rmsd_values / sum(size),
        study_ID = as.factor(study_ID)
    ) %>% 
    dplyr::group_by(x) %>% 
    dplyr::mutate(avg_rmsd_values = mean(rmsd_values, na.rm = TRUE)) %>% 
    dplyr::ungroup() %>%
    ggplot() + 
    geom_line(
        aes(x = x, y = rmsd_values), 
        linewidth = 1, alpha = 0.7, col = 'grey80'
    ) +
    geom_line(
        aes(x = x, y = avg_rmsd_values), 
        linewidth = 1, alpha = 0.7, col = 'darkred'
    ) + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 12),
        legend.position = "none"
    ) + 
    labs(
        x = "Position away from breakpoint, bp",
        y = "RMSD"
    )

ggsave(
    filename = "../figures/RMSD_pure_control.pdf",
    plot = p1,
    height = 6, width = 8
)