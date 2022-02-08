# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
k <- as.integer(args[3])
upper.limit <- as.integer(args[4])

# my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/02_Correlations/scripts/"
# breakpoint.experiment="04-Enzymatic/exp"
# k=2
# upper.limit=1
setwd(my.path)

# load dependencies
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(seqLogo))
suppressPackageStartupMessages(library(pbapply))
source("../lib/LoadKmertoneData.R")
source("../lib/GenerateKmerTable.R")
source("../lib/FilterForExtremeKmers.R")
source("../lib/LoopOverKmertoneData.R")

# --------------
# helper functions
PlotKmerClustering <- function(to.cluster, action){
  # hierarchical clustering of k-mers vs. average enrichment scores 
  if(k == 2) row.cluster.label.size=0.8
  if(k == 4) row.cluster.label.size=0.2
  if(k == 6) row.cluster.label.size=0.1
  if(k == 8) row.cluster.label.size=0.01

  colnames(to.cluster) <- paste0("exp_", seq(1, upper.limit, 1))
  df.dist <- dist(1 - to.cluster) %>% suppressWarnings()
  df.dendro <- as.dendrogram(hclust(df.dist, method = "complete"))
  transpose.df <- t(to.cluster)

  pdf(paste0("../figures/", breakpoint.experiment, "_",
            k, "-kmer_clustering_of_",
            action, "_plot.pdf"))
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
# --------------

sample.kmer <- data.table(kmer = GenerateKmerTable(k, ref = FALSE))
kmer.ref <- GenerateKmerTable(k, ref = TRUE)

# extract probability ratios from kmertone data
results <- LoopOverKmertoneData(
  upper.limit=upper.limit, 
  breakpoint.experiment=breakpoint.experiment,
  action = "ratio"
)
avg.prob <- rowMeans(results)
avg.sd <- apply(results, 1, sd)
match.ind <- match(kmer.ref$kmer, sample.kmer$kmer)

df <- data.table(
  kmer = sample.kmer$kmer[match.ind], 
  average.prob = avg.prob,
  avg.minus.sd = avg.prob-avg.sd,
  avg.plus.sd = avg.prob+avg.sd
)

setorder(df, -average.prob)
df[, kmer := forcats::fct_inorder(kmer)]

# kmer case frequency ranked bar plot
p <- df %>% 
  ggplot(aes(x = kmer, y = average.prob)) + 
  geom_bar(
    stat = "identity",
    col = "#303c54",
    fill = "#303c54"
  ) + 
  geom_errorbar(aes(
    ymin = avg.minus.sd,
    ymax = avg.plus.sd), 
    col = "darkorange"
  ) +
  geom_text(
    data = data.frame(
      xpos = Inf,
      ypos = Inf,
      annotateText = paste0(
        paste0("Top 5: ", paste0(
          FilterForExtremeKmers(df, top = TRUE, group = "average.prob"), 
        collapse = ", ")), "\n",
        paste0("Bottom 5: ", paste0(
          FilterForExtremeKmers(df, top = FALSE, group = "average.prob"), 
        collapse = ", ")), "\n",
        paste0("Most varied: ", paste0(
          FilterForExtremeKmers(df, top = TRUE, group = "avg.plus.sd"), 
        collapse = ", ")), "\n",
        paste0("Least varied: ", paste0(
          FilterForExtremeKmers(df, top = TRUE, group = "avg.minus.sd"), 
        collapse = ", "))
      ),
      hjustvar = 1.1, vjustvar = 1.1
    ),
    aes(
      x = xpos,
      y = ypos,
      hjust = hjustvar,
      vjust = vjustvar,
      label = annotateText,
      angle = 0
    ),
    size = 3
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
    title = paste0("Mean ± sd probability values of first lexicological occurring", "\n",
                   k, "-mers across all experiments"),
    x = "Kmer",
    y = "Average probability"
  )

ggsave(
  filename = paste0("../figures/", breakpoint.experiment, "_", k, "-mer_probability_vs_sd_plot.pdf"),
  plot = p,
  height = 7, width = 10
)

# average absolute z-scores vs. average probability ratios
results.zscore <- LoopOverKmertoneData(
  upper.limit=upper.limit, 
  breakpoint.experiment=breakpoint.experiment,
  action = "z-score"
)
avg.zscore <- rowMeans(results.zscore)

if(upper.limit > 1){
  avg.zscore.sd <- apply(results.zscore, 1, sd)
}

match.ind <- match(kmer.ref$kmer, sample.kmer$kmer)

df <- data.table(
  kmer = sample.kmer$kmer[match.ind], 
  average.prob = log2(avg.prob),
  average.zscore = abs(avg.zscore)
)

if(upper.limit > 1){
  df[, average.zscore.sd := avg.zscore.sd]
}

# remove any infinite values
if(any(!is.finite(rowSums(df[, 2:length(df)])))){
  df <- df[is.finite(rowSums(df[,2:length(df)])), ]  
}

if(upper.limit > 1){
  df.subset <- rbind(
    df[order(df[, average.prob], decreasing = TRUE)[1:7]],
    df[order(df[, average.prob], decreasing = FALSE)[1:7]],
    df[order(df[, average.zscore], decreasing = TRUE)[1:7]]
  )

  p <- df %>%
    ggplot() + 
    geom_point(aes(
      x = average.prob,
      y = average.zscore,
      col = ifelse(average.prob >= 0, "skyblue", "#FF7F7F")
    )) + 
    ggrepel::geom_text_repel(
      data = df.subset, 
      aes(average.prob, average.zscore, label = kmer), 
      box.padding = 0.5, 
      max.overlaps = Inf,
      min.segment.length = 0
    ) +
    scale_color_manual(
      labels = c("Depleted", "Enriched"), 
      values = c("skyblue", "#FF7F7F")
    ) + 
    coord_cartesian(
      xlim = c(-max(abs(df$average.prob))*1.1, max(abs(df$average.prob))*1.1)
    ) +
    labs(
      title = paste0(k, "-meric ratio vs. z-scores of all experiments"),
      x = "Log2 average probability ratios", 
      y = "Average absolute z-scores", 
      col = ""
    )

  ggsave(
    filename = paste0("../figures/", breakpoint.experiment, 
                      "_", k, "-kmer_ratio_vs_sd.pdf"), 
    plot = p
  )

  # zscore average vs. st.dev 
  p <- df %>%
    ggplot(aes(
      x = average.zscore, 
      y = average.zscore.sd
    )) +
    geom_point(alpha = 0.7) + 
    labs(
      title = paste0(k, "-meric average z-scores vs. standard deviations of all experiments"),
      x = "Average z-scores",
      y = "Standard deviation",
    )

  ggsave(
    filename = paste0("../figures/", breakpoint.experiment, 
                      "_", k, "-kmer_zscores_vs_sd.pdf"), 
    plot = p
  )

  for(action in c("z-score", "ratio")){
    PlotKmerClustering(to.cluster = results.zscore, action = action)
  }
}
 
if(k >= 6){
  # SeqLogos of 10% most extreme kmers
  number.of.points <- ceiling((5/100)*nrow(results))

  df.subset <- rbind(
    df[order(df[, average.prob], decreasing = TRUE)[1:number.of.points]],
    df[order(df[, average.zscore], decreasing = TRUE)[1:number.of.points]]
  )

  enriched.kmers <- c(df.subset$kmer, paste0(reverseComplement(DNAStringSet(df.subset$kmer))))

  # enriched seqlogos
  m <- consensusMatrix(enriched.kmers, as.prob = TRUE)
  p <- makePWM(m)

  pdf(paste0("../figures/", breakpoint.experiment, "_", k, "-kmer_seqlogo_info_bits_enriched.pdf")) 
  seqLogo(p, ic.scale = TRUE)
  plot.save <- dev.off()

  # enriched seqlogos
  df.subset <- rbind(
    df[order(df[, average.prob], decreasing = FALSE)[1:number.of.points]],
    df[order(df[, average.zscore], decreasing = FALSE)[1:number.of.points]]
  )

  depleted.kmers <- c(df.subset$kmer, paste0(reverseComplement(DNAStringSet(df.subset$kmer))))

  m <- consensusMatrix(depleted.kmers, as.prob = TRUE)
  p <- makePWM(m)

  pdf(paste0("../figures/", breakpoint.experiment, "_", k, "-kmer_seqlogo_info_bits_depleted.pdf")) 
  seqLogo(p, ic.scale = TRUE)
  plot.save <- dev.off()
}