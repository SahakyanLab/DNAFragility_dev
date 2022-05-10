# read arguments from job submission
# args <- commandArgs(trailingOnly = TRUE)
# my.path <- as.character(args[1])
# breakpoint.experiment <- as.character(args[2])
# k <- as.integer(args[3])
# upper.limit <- as.integer(args[4])
# action <- as.character(args[5])

my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/02_Correlations/scripts/"
k=8
action="z-score"
upper.limit=1
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

# kmer reference
sample.kmer <- data.table(kmer = GenerateKmerTable(k, ref = FALSE))
kmer.ref <- GenerateKmerTable(k, ref = TRUE)

# import kmertone data
all.files <- list.files(
  path = "../data/kmertone/", 
  pattern = paste0("score_", k, " *"), 
  recursive = TRUE
)

results <- sapply(all.files, function(file){
  out <- LoopOverKmertoneData(
    upper.limit=1,
    breakpoint.experiment=str_extract(string = file, pattern = ".+(?=/score)"),
    action=action
  )

  if(upper.limit > 1){
    out <- as.matrix(rowMeans(out))
  }

  return(out)
})

colnames(results) <- unique(str_extract(
  string = all.files, 
  pattern = "[^[:digit:]-/].+(?=/score)")
)
rownames(results) <- kmer.ref$kmer

# remove any infinite values
if(any(!is.finite(rowSums(results)))){
  results <- results[is.finite(rowSums(results)), ]
}

df.heatmap <- as_tibble(results) %>% 
  mutate(kmer = rownames(results), .before = 1) %>% 
  tidyr::gather(-kmer, key = "Variable", value = "Value") %>% 
  filter(!is.na(Value)) %>% 
  group_by(Variable) %>% 
  dplyr::mutate(Scaled = scale(Value, center = TRUE, scale = TRUE)[,1]) %>% 
  ungroup()
  
CI <- with(df.heatmap, mean(Scaled)+c(-1, 1)*1.96*sd(Scaled))

# heatmap without clustering
p <- df.heatmap %>% 
  ggplot(aes(
    x = factor(kmer), 
    y = factor(Variable), 
    fill = Scaled,
  )) + 
  geom_tile() + 
  scale_fill_gradientn(
    limits = CI,
    colors = hcl.colors(nrow(df.heatmap), "RdYlGn"),
    oob = scales::squish
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

# hierarchical clustering
df.hc <- apply(results, 2, scale, center = TRUE, scale = TRUE)
rownames(df.hc) <- rownames(results)

df.dist <- dist(1-df.hc) %>% suppressWarnings()
df.dendro <- as.dendrogram(hclust(df.dist, method = "complete"))
transpose.df <- t(df.hc)
row.clust <- as.dendrogram(hclust(dist(1-transpose.df), method = "complete"))

row.cluster.label.size <- c("2" = 0.8, "4" = 0.01, "6" = 0.01, "8" = 0.01)

cols=200

pdf(paste0("../figures/tracts_for_ml/", k, "-mer_", action, "-clustering.pdf"), height = 7, width = 12)
heatmap.2(
  transpose.df,
  offsetRow = 0,
  offsetCol = 0,
  Rowv = row.clust,
  Colv = df.dendro,
  dendrogram = "both",
  revC = FALSE,
  trace = "none",
  density.info = "none",
  col = hcl.colors(cols, "RdYlGn"),
  breaks = seq(CI[1], CI[2], length.out = cols+1),
  notecol = "black",
  cexRow = 0.4,
  cexCol = row.cluster.label.size[as.character(k)],
  labRow = rownames(transpose.df),
  labCol = colnames(transpose.df),
  margins = c(3,19),
  symkey = TRUE,
  key.xlab="Scaled z-scores"
)
plot.save <- dev.off()