# read arguments from job submission
# args <- commandArgs(trailingOnly = TRUE)
# my.path <- as.character(args[1])
# k <- as.integer(args[2])

my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/02_Correlations/scripts/"
k=4
setwd(my.path)

# load dependencies
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

# import RMSD data
all.files <- list.files(
  path = "../../01_RMSD/data/", 
  pattern = paste0("key_stats_kmer_", k, " *"), 
  recursive = TRUE 
)

results <- lapply(all.files, function(file){
  out <- fread(paste0("../../01_RMSD/data/", file))
  out <- as_tibble(out)

  if(!("curve.three" %in% colnames(out))){
    out <- cbind(out, rep(NA_real_, 4))
    colnames(out)[length(colnames(out))] <- "curve.three"
  }

  return(out)
}) 

results <- do.call(rbind, results)
colnames(results)[4:6] <- c("long.range", "short.range", "mid.range")

results <- results %>% 
  relocate(long.range, mid.range, short.range, .after = 3) %>% 
  arrange(category)

fwrite(
  results,
  paste0("../data/seq_influence/key_stats_kmer_", k, ".csv")
)

# hierarchical clustering
suppressPackageStartupMessages(library(gplots))

df <- results %>% 
  filter(rowid == "ranges") %>% 
  select(-c(category, rowid))

df.hc <- apply(df[,-1], 2, scale, center = TRUE, scale = TRUE)
rownames(df.hc) <- df$exp

df.dist <- dist(1-df.hc) %>% suppressWarnings()
df.dendro <- as.dendrogram(hclust(df.dist, method = "complete"))
transpose.df <- t(df.hc)
col.clust <- as.dendrogram(hclust(dist(1-transpose.df), method = "complete"))

cols <- 200

dir.create("../figures/seq_influence", showWarnings = FALSE)
pdf(paste0("../figures/seq_influence/", k, "-mer-clustering.pdf"), height = 10, width = 8)

# 95% confidence interval
# t.test(df.hc)
# CI <- t.test(df.hc)$conf.int

Lower.bound <- mean(df.hc, na.rm = TRUE)-2*sd(df.hc, na.rm = TRUE)
Upper.bound <- mean(df.hc, na.rm = TRUE)+2*sd(df.hc, na.rm = TRUE)

heatmap.2(
  df.hc,
  offsetRow = 0,
  offsetCol = 0,
  Rowv = df.dendro,
  Colv = FALSE,
  dendrogram = "row",
  revC = FALSE,
  trace = "none",
  density.info = "density",
  col = hcl.colors(cols, "RdYlGn"), 
  breaks = seq(Lower.bound, Upper.bound, length.out = cols+1),
  na.color = "grey",
  notecol = "black",
  cexCol = 0.6,
  cexRow = 0.6,
  labRow = rownames(df.hc),
  labCol = colnames(df.hc),
  margins = c(4,20),
  symkey = TRUE,
  key.xlab="Scaled RMSD ranges"
)
plot.save <- dev.off()