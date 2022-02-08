# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
setwd(my.path)

# load dependencies
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
source("../lib/LoadMatrixFiles.R")

# load correlation matrices
title.extract <- str_extract(string = breakpoint.experiment, pattern = "[0-9].+(?=/)")
df.ratio <- LoadMatrixFiles(action = "ratio") %>% as_tibble()
df.zscore <- LoadMatrixFiles(action = "z-score") %>% as_tibble()

colnames(df.ratio) <- str_extract(string = colnames(df.ratio), pattern = "kmer.+(?=.csv)")
colnames(df.zscore) <- str_extract(string = colnames(df.zscore), pattern = "kmer.+(?=.csv)")

# load matrix file for breakage overlap percentages
dat <- fread(
  file = paste0("../data/overlap/", title.extract, "/overlap_matrix.csv"),
  header = TRUE, 
  showProgress = FALSE
)
dat[lower.tri(dat, diag = TRUE)] <- NA
mat <- as.numeric(na.omit(unlist(dat)))

# convert data frames to long format and merge
df.ratio <- df.ratio %>% 
  mutate(df.group = "Probability ratio") %>%
  tidyr::gather(-df.group, key = "Variable", value = "Value")

df.zscore <- df.zscore %>% 
  mutate(df.group = "Z-scores") %>%
  tidyr::gather(-df.group, key = "Variable", value = "Value")

mat <- as_tibble(mat) %>% 
  mutate(df.group = "Breakage Sites") %>%
  tidyr::gather(-df.group, key = "Variable", value = "Value")

# plot distribution of correlation coefficient points
PlotDistribution <- function(d){
  ggplot() + 
    geom_density(data = d, aes(x = Value, fill = Variable), col = "black", alpha = 0.5) +
    facet_wrap(~df.group,  nrow = 3, scales = "free_x") +
    scale_x_continuous(limits = c(0.1, 1))
}

p1 <- PlotDistribution(d = mat) + 
  labs(x = "Percent Overlap", y = "") + 
  theme(legend.position = "none")

p2 <- PlotDistribution(d = df.ratio) + 
  labs(x = "", y = "") + 
    theme(legend.position = "none")

p3 <- PlotDistribution(d = df.zscore) + 
  labs(x = "Correlation Coefficient", y = "Density") +
  theme(
    legend.position = c(0.08, 0.65), 
    legend.background = element_rect(fill = NA), 
    legend.title = element_blank()
  )

ggsave(
  filename = paste0("../figures/", title.extract, "/CorrelationDistributionOverlaps.pdf"),
  plot = cowplot::plot_grid(plotlist = list(p1, p2, p3), axis = "b", ncol = 1)
)