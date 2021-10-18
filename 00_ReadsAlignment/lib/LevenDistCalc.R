suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(stringr))
pbo = pboptions(type="txt")

# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
reads.path <- as.character(args[1])
my.path <- as.character(args[2])
cores <- as.numeric(args[3])
setwd(paste0(my.path, "../lib"))

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
  
  # load files
  files <- list.files(path = paste0("../data/reads/", reads.path, "/breakpoint_positions/", chromosome, "/"),
                      pattern = "new_alignment_file")
  
  # sort files
  files <- str_sort(files, numeric = TRUE)
  
  tables <- lapply(files, function(x){
    fread(paste0("../data/reads/", reads.path, "/breakpoint_positions/", chromosome, "/", x), 
               sep = ",", header = TRUE)
  })
  
  # concatenate files and change class to numeric
  df <- do.call(rbind , tables) %>% 
    as_tibble() %>%
    rename(start.pos = bp.start.pos,
           lev.dist  = lev.dist) %>%
    drop_na()
  
  # convert any non-integer columns to integers
  df <- df %>% 
    mutate(across(where(is.character), as.integer)) %>%
    suppressWarnings()
  
  # remove any NAs
  df <- df %>%
    drop_na()
  
  # obtain levenshtein frequency table
  lev.dist.df <- df %>%
    select(-start.pos) %>%
    group_by(lev.dist) %>%
    summarise(count = sum(freq))
  
  # return frequency average levenshtein distance
  Mean <- sum(lev.dist.df$lev.dist*lev.dist.df$count)/sum(lev.dist.df$count)
  SD <- sqrt(sum((lev.dist.df$lev.dist)**2*lev.dist.df$count)/(sum(lev.dist.df$count)-1))
  
  return(list(Mean, SD))
}, cl = cores)

# concatenate results and calculate overall statistics
df <- do.call(rbind, lev.dist) %>%
  as_tibble() %>%
  mutate(across(where(is.list), as.numeric)) %>%
  rename(Mean = V1, SD = V2) %>%
  mutate(SD1 = Mean+SD,
         SD2 = Mean+2*SD)

# rename column
names <- paste0("Chr", sort(seq(1:22), decreasing = TRUE))

# obtain bar plot of levenshtein distance
lev.plot <- df %>%
  rownames_to_column() %>%
  mutate(across(where(is.character), as.numeric)) %>%
  arrange(desc(rowname)) %>%
  mutate(rowname = names) %>%
  rename(Chromosomes = rowname) %>%
  mutate(Chromosomes = forcats::fct_inorder(Chromosomes)) %>%
  ggplot() +
  geom_bar(aes(x = Chromosomes,
               y = Mean),
           stat = "identity",
           fill = "skyblue",
           alpha = 0.7) + 
  geom_pointrange(aes(x = Chromosomes,
                      y = Mean,
                      ymin = Mean,
                      ymax = SD1),
                  color = "orange",
                  alpha = 1,
                  size = 0.8) +
  coord_flip() +
  labs(title = "Average Levenshtein Distance plotted as Mean + 1 St.Dev", 
       subtitle = paste0("Overall Average: ", round(mean(df$Mean), 3)))

ggsave(paste0("../figures/", reads.path, "/AvgLevenshteinDistance.pdf"),
       plot = lev.plot)

# save lev dist 
fwrite(x = df, 
       file = paste0("../data/average_levdist/", reads.path, "/AvgLevenshteinDistance.csv"), 
       row.names = FALSE)