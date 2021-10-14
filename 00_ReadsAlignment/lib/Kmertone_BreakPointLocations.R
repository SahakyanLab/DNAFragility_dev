suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(data.table))

# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
reads_path <- as.character(args[1])
my_path <- as.character(args[2])

setwd(paste0(my_path, "../lib/"))

# ---------------
pb <- txtProgressBar(min = 1, max = 23, style = 3)

# obtain average levenshtein distance for all chromosomes
for(i in 1:22){
  if(i == 23){
    chromosome = "chrX"
  } else if (i == 24) {
    chromosome = "chrY"
  } else {
    chromosome = paste0("chr", i)
  }
  
  # load data sets
  files <- list.files(path = paste0("../../data/reads/", reads_path, "/breakpoint_positions/", chromosome, "/"),
                      pattern = "alignment_file_")
  
  # concatenate files
  tables <- lapply(files, function(x){
    fread(file = paste0("../../data/reads/", reads_path, "/breakpoint_positions/",
                        chromosome,"/", x), sep = ",", header = TRUE)
  })
  df <- do.call(rbind, tables) %>% 
    as_tibble()
  
  # convert any non-integer columns to integers
  df <- df %>% 
    mutate(across(where(is.character), as.integer)) %>%
    suppressWarnings()
  
  # remove any NAs
  df <- df %>%
    drop_na()
  
  # save bottom ~95% of the levenshtein distance score
  # import average lev dist
  lev.dist.df <- read.table(file = paste0("../../data/two_mers/", reads_path, "/AvgLevenshteinDistance.csv"),
                            sep = ",",
                            header = TRUE) %>% as_tibble()
  
  df <- df %>%
    filter(lev_dist < mean(lev.dist.df$SD1))
  
  # keep only bp start position column
  chr.name = rep(chromosome, dim(df)[1])
  df <- df %>% 
    select(bp_start_pos) %>% 
    mutate("chromosome" = chr.name)
  
  # save start breakpoint positions and chromosome name as txt file
  write.table(x = df, row.names = FALSE, file = paste0("../../data/kmertone/breakpoints/", chromosome, ".txt"))

  setTxtProgressBar(pb, i)
}