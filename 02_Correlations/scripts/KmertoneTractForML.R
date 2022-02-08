# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
k <- as.integer(args[3])
upper.limit <- as.integer(args[4])
action <- as.character(args[5])
setwd(my.path)

# load dependencies
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(pbapply))
source("../lib/GenerateKmerTable.R")
source("../lib/LoadKmertoneData.R")
source("../lib/LoopOverKmertoneData.R")

# ----------------
# helper functions
PlotKmerClustering <- function(to.cluster, action){
  if(k == 2) row.cluster.label.size=0.8
  if(k == 4) row.cluster.label.size=0.2
  if(k == 6) row.cluster.label.size=0.1
  if(k == 8) row.cluster.label.size=0.01

  df.dist <- dist(1 - to.cluster) %>% suppressWarnings()
  df.dendro <- as.dendrogram(hclust(df.dist, method = "complete"))
  transpose.df <- t(to.cluster)

  pdf(paste0("../figures/tracts_for_ml/", k, "-mer_", action, "-clustering.pdf"))
  heatmap.2(
    transpose.df,
    Rowv = FALSE,
    Colv = df.dendro,
    dendrogram = "col",
    revC = FALSE,
    trace = "none",
    density.info = "none",
    notecol = "black",
    cexRow = 0.8,
    cexCol = row.cluster.label.size,
    labRow = rownames(transpose.df),
    labCol = colnames(transpose.df)
  )
  plot.save <- dev.off()
}
# ----------------

# kmer reference
sample.kmer <- data.table(kmer = GenerateKmerTable(k, ref = FALSE))
kmer.ref <- GenerateKmerTable(k, ref = TRUE)

# import kmertone data
all.files <- list.files(
  path = "../data/kmertone/", 
  pattern = paste0("score_", k, " *"), 
  recursive = TRUE
)

df <- tibble(
  breakpoint.experiment = str_extract(
    string = all.files, 
    pattern = ".+(?=_[:digit:]+/)")) %>% 
  group_by(breakpoint.experiment) %>% 
  summarise(ind = n())

results <- sapply(1:nrow(df), function(ind){
  upper.limit <- pull(df[ind, "ind"])

  out <- LoopOverKmertoneData(
    upper.limit=upper.limit,
    breakpoint.experiment=pull(df[ind, "breakpoint.experiment"]),
    action=action
  )

  if(upper.limit > 1){
    out <- as.matrix(rowMeans(out))
  }

  return(out)
})

colnames(results) <- unique(str_extract(string = all.files, pattern = "[^[:digit:]-/]+"))
rownames(results) <- kmer.ref$kmer

# remove any infinite values
if(any(!is.finite(rowSums(results)))){
  results <- results[is.finite(rowSums(results)), ]
}

# hiarchical clustering
PlotKmerClustering(to.cluster = results, action = action)

# heatmap without clustering
p <- as_tibble(results) %>% 
  mutate(kmer = rownames(results), .before = 1) %>% 
  tidyr::gather(-kmer, key = "Variable", value = "Value") %>% 
  ggplot(aes(
    x = kmer, 
    y = Variable, 
    fill = Value
  )) + 
  geom_tile() +
  scale_fill_gradientn(
    colors = hcl.colors(nrow(results), "YlOrRd")
  ) + 
  theme_bw() + 
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + 
  labs(
    x = "Kmers",
    y = ""
  )

ggsave(
  filename = paste0("../figures/tracts_for_ml/", k, "-mer_", action, "-nocluster.pdf"),
  plot = p,
  height = 7, width = 12
)