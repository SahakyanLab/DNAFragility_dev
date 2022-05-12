# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
k <- as.integer(args[2])
auto.fit <- as.logical(args[3])
setwd(my.path)

# load dependencies
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

# import RMSD data
all.files <- list.files(
  path = "../../01_RMSD/data/", 
  pattern = paste0(ifelse(auto.fit, "key", "new"), 
                   "_stats_kmer_", k, " *"), 
  recursive = TRUE 
)

if(auto.fit) curve.labels <- c("long.range", "mid.range", "short.range")

results <- lapply(all.files, function(file){
  out <- fread(paste0("../../01_RMSD/data/", file))
  out <- as_tibble(out)

  if(auto.fit){
    if(!("curve.three" %in% colnames(out))){
      out <- cbind(out, rep(NA_real_, 4))
      colnames(out)[length(colnames(out))] <- "curve.three"
    }
  } else {
    missing <- which(!curve.labels %in% colnames(out))
    if(any(missing)){
      out <- cbind(out, rep(NA_real_, 4))
      colnames(out)[length(colnames(out))] <- curve.labels[missing]
    }
  }

  return(out)
}) 

results <- do.call(rbind, results)
if(auto.fit) colnames(results)[4:6] <- c("long.range", "short.range", "mid.range")

results <- results %>% 
  relocate(long.range, mid.range, short.range, .after = 3) %>% 
  arrange(category)

fwrite(
  results,
  paste0("../data/seq_influence/", 
         ifelse(auto.fit, "key", "new"), 
         "_stats_kmer_", k, ".csv")
)

# hierarchical clustering
suppressPackageStartupMessages(library(gplots))

df <- results %>% 
  filter(rowid == "ranges") %>% 
  select(-c(category, rowid))

if(!auto.fit){
  df <- df %>% 
    dplyr::relocate(short.range, mid.range, long.range, .after = 1)
}

df.hc <- apply(df[,-1], 2, scale, center = FALSE, scale = FALSE)
rownames(df.hc) <- df$exp

df.dist <- dist(1-df.hc) %>% suppressWarnings()
df.dendro <- as.dendrogram(hclust(df.dist, method = "complete"))
transpose.df <- t(df.hc)
col.clust <- as.dendrogram(hclust(dist(1-transpose.df), method = "complete"))

dir.create("../figures/seq_influence", showWarnings = FALSE)
pdf(paste0("../figures/seq_influence/", k, "-mer-clustering", 
           ifelse(auto.fit, "auto_fit", "optimised"),
           ".pdf"), height = 10, width = 8)

cols <- 200
if(auto.fit){
  heatmap.breaks <- seq(min(df.hc, na.rm = TRUE), 
                        max(df.hc, na.rm = TRUE), 
                        length.out = cols+1)
} else {
  Lower.bound <- mean(df.hc, na.rm = TRUE)-2*sd(df.hc, na.rm = TRUE)
  Upper.bound <- mean(df.hc, na.rm = TRUE)+2*sd(df.hc, na.rm = TRUE)

  heatmap.breaks <- seq(Lower.bound, 
                        Upper.bound, 
                        length.out = cols+1)
}

heatmap.2(
  df.hc,
  offsetRow = 0,
  offsetCol = 0,
  Rowv = df.dendro,
  Colv = FALSE,
  dendrogram = "row",
  revC = FALSE,
  trace = "none",
  density.info = "histogram",
  col = hcl.colors(cols, "RdYlGn"), 
  breaks = heatmap.breaks,
  na.color = "grey",
  notecol = "black",
  cexCol = 0.6,
  cexRow = 0.6,
  labRow = rownames(df.hc),
  labCol = colnames(df.hc),
  margins = c(4,20),
  key.xlab = "RMSD ranges"
)
plot.save <- dev.off()

if(!auto.fit){
  # Plot individual distributions for each range
  p <- results %>% 
      filter(rowid == "ranges") %>% 
      select(exp, category, contains("range")) %>% 
      tidyr::gather(-c(exp,category), key = "Ranges", value = "Value") %>% 
      mutate(
        Ranges = factor(Ranges, 
        levels = c("short.range", "mid.range", "long.range"))
      ) %>% 
      ggplot(aes(x = Value)) + 
      geom_histogram(aes(
        fill = category),
        bins = 50,
        colour = "black"
      ) +
      facet_wrap(~Ranges, scales = "free_x") +
      labs(
        title = "Sequence influences clustered into 3 separate ranges",
        subtitle = paste0("Kmer: ", k),
        x = "Ranges",
        y = "Count"
      )

  suppressWarnings(ggsave(
    plot = p, 
    filename = paste0("../figures/seq_influence/", k, "-mer-ranges.pdf"),
    height = 7, width = 13
  ))
}