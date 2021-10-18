suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(stringr))
pbo = pboptions(type="txt")

# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
reads.path <- as.character(args[1])
my.path <- as.character(args[2])
cores <- as.numeric(args[3])

setwd(paste0(my.path, "../lib/"))

# ---------------
# obtain average levenshtein distance for all chromosomes
lev.dist <- pblapply(1:22, function(i){
  if(i == 23){
    chromosome = "chrX"
  } else if (i == 24) {
    chromosome = "chrY"
  } else {
    chromosome = paste0("chr", i)
  }
  
  # load data sets
  files <- list.files(path = paste0("../data/reads/", reads.path, "/breakpoint_positions/", chromosome, "/"),
                      pattern = "new_alignment_file")
  
  # sort files
  files <- str_sort(files, numeric = TRUE)
  
  # concatenate files
  tables <- lapply(files, function(x){
    fread(file = paste0("../data/reads/", reads.path, "/breakpoint_positions/",
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
  lev.dist.df <- fread(file = paste0("../data/average_levdist/", reads.path, "/AvgLevenshteinDistance.csv"),
                       sep = ",",
                       header = TRUE) %>% as_tibble()
  
  df <- df %>%
    filter(lev.dist < mean(lev.dist.df$SD1))
  
  # keep only bp start position column
  chr.name = rep(chromosome, dim(df)[1])
  df <- df %>% 
    select(bp.start.pos) %>% 
    mutate("chromosome" = chr.name)
  
  # save start breakpoint positions and chromosome name as txt file
  fwrite(x = df, row.names = FALSE, 
         file = paste0("../../data/kmertone/", reads.path, "/breakpoints/", chromosome, ".txt"))
}, cl = cores)