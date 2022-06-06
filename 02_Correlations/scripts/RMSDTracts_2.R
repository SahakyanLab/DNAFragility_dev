# read arguments from job submission
# args <- commandArgs(trailingOnly = TRUE)
# my.path <- as.character(args[1])
# k <- as.integer(args[2])

# By range, not by kmer
my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/02_Correlations/scripts/"
setwd(my.path)

# load dependencies
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

# import RMSD data
all.files <- list.files(
  path = "../../01_RMSD/data/", 
  pattern = "new_stats_kmer_*", 
  recursive = TRUE 
)

files.kmers <- as.integer(str_extract(
    string = basename(all.files), 
    pattern = "[:digit:]"
))
curve.labels <- c("long.range", "mid.range", "short.range")

results <- lapply(1:length(all.files), function(i){
  out <- fread(paste0("../../01_RMSD/data/", all.files[i]))
  out <- as_tibble(out) %>% mutate(kmer = files.kmers[i])

  missing <- which(!curve.labels %in% colnames(out))
  if(any(missing)){
    out <- cbind(out, rep(NA_real_, 4))
    colnames(out)[length(colnames(out))] <- curve.labels[missing]
  }

  return(out)
}) 

results <- do.call(rbind, results)

results <- results %>% 
  filter(rowid == "contribution") %>% 
  relocate(long.range, mid.range, short.range, .after = 3) %>% 
  arrange(category)

categories <- c("Mechanical", "Biological", "Cell-free", "Ancient", "Enzymatic")

k=4
results.long <- lapply(categories, function(x){
  return(
    results %>% 
      filter(kmer == k) %>% 
      select(-c(exp, rowid, kmer)) %>% 
      filter(category == x) %>% 
      tidyr::gather(-category, key="range", value="value") %>% 
      mutate(category = x)
  )
})
results.long <- do.call(rbind, results.long)

p <- results.long %>% 
  ggplot(aes(
    x = value,
    y = ..density..
  )) + 
  geom_histogram(
    aes(x = value),
    bins = 30,
    fill = "skyblue",
    col = "black"
  ) + 
  facet_grid(
    vars(range), 
    vars(category)
  ) + 
  labs(
    title = paste0("Kmer: ", k),
    x = "Percent Contribution",
    y = "Density"
  )

ggsave(
  plot = p,
  filename = paste0("../figures/Density-Percent-Contribution_kmer_", k, ".pdf"),
  height = 7, width = 10
)