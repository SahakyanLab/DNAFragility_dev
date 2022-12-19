# read arguments from job submission
# args <- commandArgs(trailingOnly = TRUE)
# my.path <- as.character(args[1])
# ind <- as.numeric(args[2])
# control <- as.logical(args[3])
# seed <- as.numeric(args[4])
# k <- as.numeric(args[5])

# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(cluster)))
suppressPackageStartupMessages(suppressWarnings(library(doParallel)))
suppressPackageStartupMessages(suppressWarnings(library(foreach)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults

my.path="/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints_v2/03_FitCurves/scripts"
setwd(my.path)
source("../lib/SequenceEffect.R")

df <- fread("../../data/org_file.csv")
df <- df[`DSB Map` == "TRUE"]
ind=26
exp.name <- paste0(
    df[ind, `Fragmentation type`], "/",
    df[ind, `Experiment folder`]
)
seed=1234
k=8
from_file=TRUE
cores=ceiling(parallel::detectCores()/2)

# begin preprocessing steps
set.seed(seed)
rng <- sample(100000, size = 10, replace = FALSE)
for(seeds in rng){
    seq_effect <- SequenceEffect$new(
        chr = 1,
        which_exp_ind = ind,
        control = TRUE,
        seed = seeds,
        cores = cores
    )
    seq_effect$calc_seq_effect(
        k = k, 
        rmsd.range = c(-301, 301),
        from_file = from_file
    )
    seq_effect$quick_plot_check()
    print(seq_effect$plots[[paste0("kmer_", k)]])
}

# statistical test of significance
true_bp <- readRDS(paste0("../data/", exp.name, "/rmsd_kmer_", k,".Rdata"))
false_bp.files <- list.files(
    path = paste0("../data/", exp.name),
    pattern = paste0("^control_rmsd_kmer_", k),
    full.names = TRUE
)
false_bp <- lapply(false_bp.files, function(x){
    readRDS(x)
})

# statistical test for each point
results <- lapply(1:length(true_bp), function(x){
    false_bp_first_val <- sapply(false_bp, `[[`, x)
    stat.test <- wilcox.test(
        false_bp_first_val, 
        mu = true_bp[x],
        conf.level = 0.95
    )
    stat.test.pvalue <- stat.test$p.value
    stat.signif <- stat.test.pvalue < 0.05
    return(list(stat.test.pvalue, stat.signif))
})
stat.test.pvalue <- sapply(results, `[[`, 1)
stat.signif <- sapply(results, `[[`, 2)
all(stat.signif)

# plot results of control breakpoints
false_bp_mean <- sapply(1:length(true_bp), function(x){
    mean(sapply(false_bp, `[[`, x), na.rm = TRUE)
})
false_bp_df <- lapply(1:length(false_bp.files), function(x){
    temp <- readRDS(false_bp.files[x])
    limits <- length(temp)/2-1
    temp <- data.table(
        x = -limits:(length(temp)-limits-1),
        y = temp,
        exp = as.factor(paste0("exp_", x)),
        bp.class = FALSE
    )
    return(temp)
})
false_bp_df <- rbindlist(false_bp_df)

# true breakpoints plot
limits <- length(true_bp)/2-1
true_bp_df <- tibble(
    x = -limits:(length(true_bp)-limits-1),
    y = true_bp,
    exp = "exp_true",
    bp.class = TRUE)

# false breakpoints mean
limits <- length(false_bp_mean)/2-1
false_bp_mean_df <- tibble(
    x = -limits:(length(false_bp_mean)-limits-1),
    y = false_bp_mean,
    exp = "exp_false",
    bp.class = FALSE)

bp_df_all <- rbind(false_bp_df, true_bp_df)

p1 <- as_tibble(bp_df_all) %>%
    ggplot(aes(x = x, y = y)) +
    geom_line(
        linewidth = 0.8,
        alpha = 0.3) +
    geom_line(
        data = false_bp_mean_df,
        aes(x = x, y = y),
        linewidth = 0.8) + 
    geom_line(
        data = true_bp_df,
        aes(x = x, y = y),
        linewidth = 0.8) + 
    facet_grid(~bp.class, scales = "free_y") + 
    theme_bw() + 
    labs(
        x = "Position away from breakpoint",
        y = "RMSD",
        title = paste0("Exp: ", exp.name),
        subtitle = paste0(
            "Kmer: ", k, ".\n",
            "All statistically significant from true breakpoints? ", 
            all(stat.signif),  ".\n",
            "Seeds used: ",
            paste0(rng, collapse = " ", sep=""),
            ".")
    )

dir.create(
    path = paste0("../figures/", exp.name, "/"), 
    showWarnings = FALSE,
    recursive = TRUE
)
ggsave(
    filename = paste0("../figures/", exp.name, "/", 
                      "true_vs_control_bp_kmer_", k, 
                      ifelse(from_file, "_WITH-SEQ-INFLUENCE.pdf", 
                      ".pdf")),
    plot = p1,
    height = 8,
    width = 10
)