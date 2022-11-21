my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/01_RMSD/scripts/"

# breakpoint.experiment="21-BLISS/CD34_Top2_mediated_DSBs/ETO_rep1"
# category="Biological"

breakpoint.experiment="00-Ultrasonication/Simons_exp_1"
category="Mechanical"

setwd(my.path)

suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))

source("../lib/FitGMM.R")
source("../lib/IntegralFunctions.R")
source("../lib/CalcCI.R")
source("../lib/MakePlot.R")

# load data sets
files <- list.files(
  path = paste0("../data/", breakpoint.experiment, "/all_chr/"),
  pattern = "-rmsd_kmer_"
)
files <- str_sort(files, numeric = TRUE)
# extract chromosome number
chr.nr <- str_extract(string = files, pattern = "chr_[:digit:]*")
# extract kmer 
kmer.nr <- str_extract(string = files, pattern = "kmer_[:digit:]*")

data.sets <- lapply(1:length(files), function(i){
  out <- readRDS(
    file = paste0("../data/", breakpoint.experiment, "/all_chr/", files[i])
  )

  limits <- length(out)/2-1
  out <- as_tibble(out) %>% 
    mutate(
      x = -limits:(length(out)-limits-1),
      kmer = as.factor(kmer.nr[i]),
      chr = as.factor(chr.nr[i])
    ) %>% 
    dplyr::rename(y = value)

  return(out)
})
data.sets <- do.call(rbind, data.sets)

# Overall RMSD plots
p <- data.sets %>%
    # filter(kmer != "kmer_10") %>% 
    ggplot(aes(
      x = x, 
      y = y
    )) + 
    geom_line(size = 0.6) + 
    facet_grid(
      vars(kmer),
      vars(chr),
      scales = "free_y"
    ) + 
    theme_bw() + 
    labs(
        title = "RMSD profiles for 22 chromosomes",
        x = "Position away from breakpoint",
        y = "RMSD"
    )

ggsave(
    filename = paste0("../figures/", breakpoint.experiment, 
                        "/all_chr-sidebyside_RMSD.pdf"),
    plot = p, 
    height = 10, width = 19
)