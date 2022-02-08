# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
experiment.num <- as.character(args[3])
chromosome <- as.character(args[4])
setwd(my.path)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
source("../lib/CalcCI.R")
source("../lib/FitGMM.R")
source("../lib/MakePlot.R")

# load data sets
files <- list.files(
  path = paste0("../data/", breakpoint.experiment, "_", experiment.num),
  pattern = "rmsd_kmer_"
)
file.name <- paste0("kmer_", seq(from = 2, to = 10, by = 2))
files <- str_sort(files, numeric = TRUE)

data.sets <- lapply(1:length(files), function(x){
  out <- readRDS(
    file = paste0("../data/", breakpoint.experiment, "_", experiment.num, "/", files[x])
  )

  out <- as_tibble(out) %>% 
    mutate(kmer = as.factor(file.name[x]), 
           x = -300.5:(length(out)-300.5-1)) %>% 
    rename(y = value)
  return(out)
})

data.sets <- do.call(rbind, data.sets)

# Overall RMSD plots
plots <- data.sets %>%
  ggplot(aes(x = x, y = y)) + 
  geom_line(size = 0.8) + 
  facet_wrap(~kmer, ncol = 2, scales = "free_y") + 
  theme_bw() + 
  labs(
    x = "Position away from breakpoint",
    y = "RMSD"
  )

ggsave(
  filename = paste0("../figures/", breakpoint.experiment, "_", experiment.num,
                    "/chr", chromosome, "_sidebyside_RMSD.pdf"),
  plot = plots, 
  height = 7, width = 7
)

# GMM plots
data.sets <- data.sets %>% 
  mutate(across(where(is.factor), as.character)) %>% 
  filter(kmer != "kmer_2" & kmer != "kmer_4")

p <- vector(mode = "list", length = 3)
unique.k <- unique(data.sets$kmer)

for(k in 1:length(unique.k)){
  dat <- data.sets %>% filter(kmer == unique.k[k])
  curvefits <- FitGMM(dat = dat, ind = unique.k[k], sigma3 = 3)
  p[[k]] <- MakePlot(dat = dat, k = unique.k[k], curve.vals = curvefits)
  # if(k == 1) print(p) # temp. avoid an unsolved ggplot bug
}

ggsave(
  filename = paste0("../figures/", breakpoint.experiment, "_", experiment.num,
                    "/chr", chromosome, "_RMSD_all_GMMFit.pdf"),
  plot = cowplot::plot_grid(plotlist = p, axis = "b", ncol = 1), 
  height = 15, width = 7
)
