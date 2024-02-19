suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(gridExtra)))
suppressPackageStartupMessages(suppressWarnings(library(cluster)))
suppressPackageStartupMessages(suppressWarnings(library(mclust)))
suppressPackageStartupMessages(suppressWarnings(library(gplots)))
data.table::setDTthreads(threads = 1)

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

# vars
kmer=8
action="zscore"

# import processed kmertone data
data.set <- fread(
  paste0(
    "../data/kmertone/QueryTable/Extended_QueryTable_kmer-",
    kmer, "_", action, ".csv"
  ),
  showProgress = FALSE
)

org_file <- fread("../../data/org_file.csv", showProgress = FALSE)
org_file[, exp := paste0(`Fragmentation type`, "/", `Experiment folder`)]

# label categories
inds <- match(data.set$category, org_file$Category_Main)
labels <- org_file$Category_general[inds]
labels[which(is.na(labels))] <- "QM_parameters"
gen.class <- org_file$Class[inds]
gen.class[which(is.na(gen.class))] <- "QM_parameters"
category.hex <- org_file$Category_Hex[inds]

rows.to.keep <- which(grepl(pattern = "DNA_Breakage", gen.class))
data.set <- data.set[rows.to.keep,]
gen.class <- gen.class[rows.to.keep]
category.labels <- labels[rows.to.keep]
category.hex <- category.hex[rows.to.keep]
data.set.labels <- data.set$category

data.set.kmers <- colnames(data.set[,-"category"])
df.norm <- apply(data.set[,-"category"], 1, scale, 
                 scale = TRUE, center = TRUE)
df.norm <- t(df.norm)
rownames(df.norm) <- data.set.labels
colnames(df.norm) <- data.set.kmers

# remove any infinite values
df.norm <- df.norm[,is.finite(colSums(df.norm))]

########################################################################
# Explore PCA
########################################################################
# compute the PCA of Data
my_pca <- prcomp(
  t(df.norm), 
  retx = TRUE, 
  center = FALSE, 
  scale. = FALSE
)

# Cumulative and proportion of variance explained (CVE, PVE)
variance <- my_pca$sdev^2 
prop_variance <- variance/sum(variance)
cumulative_variance <- cumsum(prop_variance)

pc_cumvar <- tibble(
  PC = 1:length(variance),
  PVE = prop_variance,
  CVE = cumulative_variance
)
  
p1 <- pc_cumvar %>%  
  tidyr::pivot_longer(
    cols = !PC,
    names_to = "key",
    values_to = "value"
  ) %>% 
  ggplot(aes(x = PC, y = value)) + 
  geom_line() + 
  geom_point(size = 2) + 
  theme_bw() + 
  theme_classic() + 
  theme(text = element_text(size = 20)) + 
  facet_wrap(~key, ncol = 1, scales = "free_y") + 
  labs(y = "Variance Explained")

dir.create("../figures/PCA/", showWarnings = FALSE)
ggsave(
    filename = paste0("../figures/PCA/",
                      "CVE_kmer-", kmer, "_", 
                      action, ".pdf"), 
    plot = p1,
    width = 8,
    height = 8
)

# PC1 and PC2 feature loadings
feature.loading <- as_tibble(my_pca$rotation[, 1:2]) %>%
  dplyr::mutate(
    feature = rownames(my_pca$rotation[, 1:2]),
    .before = 1) %>% 
  tidyr::pivot_longer(
    cols = !feature,
    names_to = "key", 
    values_to = "value"
  ) %>% 
  dplyr::group_by(key) %>% 
  dplyr::arrange(value, .by_group = TRUE)

loading.key <- stringr::str_sort(unique(feature.loading$key), numeric = TRUE)
p1 <- lapply(loading.key, function(pca){
  label <- switch(pca,
    "PC1" = paste0(
      "Feature Loading PC1 (", 
      signif((pc_cumvar$PVE*100)[1], 3), "%)"
    ),
    "PC2" = paste0(
      "Feature Loading PC2 (", 
      signif((pc_cumvar$PVE*100)[2], 3), "%)"
    )
  )

  p1 <- feature.loading %>% 
    dplyr::filter(key == pca) %>% 
    dplyr::slice_tail(n = 60) %>% 
    dplyr::mutate(feature = forcats::fct_inorder(feature)) %>% 
    ggplot(aes(x = value, y = feature)) + 
    geom_segment(
      aes(xend = min(value), yend = feature), 
      col = "grey80",
      linewidth = 1.1
    ) + 
    geom_point(col = "orange", size = 2) + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    labs(x = "", y = "", title = label)
  return(p1)
})

png(
  paste0("../figures/PCA/Feature_loadings_PC1_PC2", 
          "_kmer-", kmer, "_", action, ".png"),
    height = 10, width = 16, 
    units = "in", res = 300
)
do.call(gridExtra::grid.arrange, c(p1, ncol = 2))
plots.saved <- dev.off()

# extract PC axes for plotting
PCAloadings <- tibble(
  Variables = rownames(my_pca$rotation),
  Category = category.labels,
  PC1 = my_pca$rotation[, "PC1"],
  PC2 = my_pca$rotation[, "PC2"],
  cols = category.hex
)

# without labels
p1 <- PCAloadings %>% 
  dplyr::mutate(
    Variables = stringr::str_replace_all(
      string = Variables,
      pattern = "_",
      replacement = " ")) %>% 
    # Variables = stringr::str_replace(
    #   string = Variables,
    #   pattern = "^(\\S+) (\\S+)",
    #   replacement = "\\1 \\2\n")) %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(
    x = 0, y = 0, 
    xend = PC1, yend = PC2),
    arrow = arrow(length = unit(1/2, "picas")),
    col = category.hex) +
  ggrepel::geom_text_repel( 
    aes(x = PC1, y = PC2, label = Variables),
    col = PCAloadings$cols,
    box.padding = 0.5, 
    max.overlaps = Inf,
    min.segment.length = unit(0.1, "lines"),
    segment.size = 0.2,
    segment.color = PCAloadings$cols,
    size = 5.5
  ) +
  theme_bw() + 
  theme_classic() + 
  theme(text = element_text(size = 25)) +
  labs(
    x = paste0("PC1 (", signif(pc_cumvar$PVE[1]*100, 3)," %)"),
    y = paste0("PC2 (", signif(pc_cumvar$PVE[2]*100, 3)," %)")
  )

ggsave(
    filename = paste0("../figures/PCA/",
                      "BreakageType_AllLabels_kmer-", kmer, "_", 
                      action, ".pdf"), 
    plot = p1,
    width = 18,
    height = 15
)

# with labels
PCAloadings.uniquelabels <- PCAloadings %>% 
  dplyr::filter(
    Variables == "AsiSI_restriction_enzyme_AID-Diva_cells" |
    Variables == "AsiSI_restriction_enzyme_K562_cells" |
    Variables == "AsiSI_restriction_enzyme_MCF10A_cells" |
    Variables == "EcoRV_restriction_enzyme_HeLa_cells" |
    Variables == "Nt_BbvCI_restriction_enzyme_K562_cells" |
    Variables == "Twist_library_C25cl48_cells" |
    Variables == "Cancer_Esophageal" |
    Variables == "Healthy" |
    Variables == "Autoimmune_Systemic_Lupus" |
    Variables == "HCT116_cells_shWRN_WRN_loss"
  ) %>% 
  dplyr::group_by(Category) %>% 
  dplyr::mutate(count = 1:(dplyr::n()), .after = Category) %>% 
  dplyr::mutate(
    key.label = ifelse(
      Category == "Enzymatic", paste0("E", count),
      ifelse(Category == "Biological", paste0("B", count), 
      paste0("CF", count)))
  ) %>% 
  dplyr::ungroup()

p2 <- PCAloadings %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(
    x = 0, y = 0, 
    xend = PC1, yend = PC2),
    arrow = arrow(length = unit(1/2, "picas")),
    col = category.hex,
    alpha = 0.7) +
  ggrepel::geom_text_repel(
    data = PCAloadings.uniquelabels, 
    aes(x = PC1, y = PC2, label = key.label),
    col = PCAloadings.uniquelabels$cols,
    box.padding = 0.5, 
    max.overlaps = Inf,
    min.segment.length = unit(0.1, "lines"),
    segment.size = 0.2,
    segment.color = PCAloadings.uniquelabels$cols,
    size = 7
  ) +
  theme_bw() + 
  theme_classic() + 
  theme(text = element_text(size = 20)) +
  labs(
    x = paste0("PC1 (", signif(pc_cumvar$PVE[1]*100, 3)," %)"),
    y = paste0("PC2 (", signif(pc_cumvar$PVE[2]*100, 3)," %)")
  )

ggsave(
    filename = paste0("../figures/PCA/",
                      "BreakageType_SimpleLabels_kmer-", 
                      kmer, "_", action, ".pdf"), 
    plot = p2,
    width = 9,
    height = 7
)

# # clustering PCA results
# feat.load.wide <- feature.loading %>% 
#   tidyr::spread(key, value) %>% 
#   tibble::column_to_rownames("feature") %>%
#   ggrepel::geom_text_repel(
#     data = PCAloadings.uniquelabels, 
#     mapping = aes(x = PC1, y = PC2, label = key.label),
#     col = PCAloadings.uniquelabels$cols,
#     box.padding = 0.5, 
#     max.overlaps = Inf,
#     min.segment.length = unit(0.1, "lines"),
#     segment.size = 0.2,
#     segment.color = PCAloadings.uniquelabels$cols,
#     size = 7
#   )
# factoextra::eclust(feat.load.wide, "kmeans", hc_metric = "euclidean", k = 6)

################################################################################
# Cluster PC1 and PC2 by breakage type
my_pca.clust.breaktype <- data.table(
  breakage_category = rownames(my_pca$rotation[, 1:2]),
  PC1 = my_pca$rotation[, 1],
  PC2 = my_pca$rotation[, 2]
)

# calculate within groups sum of squares
wss <- (nrow(my_pca.clust.breaktype)-1)*sum(apply(
  my_pca.clust.breaktype[, -"breakage_category"],2,var
))

set.seed(1234)
for(i in 2:15){
  wss[i] <- sum(kmeans(
    my_pca.clust.breaktype[, -"breakage_category"], 
    centers = i)$withinss
  )
}

# finding optimal nr of clusters based on BIC
set.seed(1234)
d_clust <- Mclust(
  my_pca.clust.breaktype[, -"breakage_category"], 
  G=1:15, 
  verbose = FALSE
)
m.best <- dim(d_clust$z)[2]

# plot combined results
pdf(paste0("../figures/PCA/",
           "Clustering_PCs_BreakageType_kmer-", 
           kmer, "_", action, ".pdf"),
  height = 7, width = 12
)
par(mfrow = c(1,2))
plot.Mclust(d_clust, what = "classification", main = TRUE)
plot(1:15, wss, 
  type="b", 
  xlab="Number of Clusters",
  ylab="Within groups sum of squares"
)
plots.saved <- dev.off()

# label clusters
set.seed(1234)
custom_clust <- Mclust(
  my_pca.clust.breaktype[, -"breakage_category"], 
  G=6, 
  verbose = FALSE
)
my_pca.clust.breaktype[, kmeans.clust := custom_clust$classification]

pdf(paste0("../figures/PCA/",
           "Clustering_PCs_BreakageType_MY_MCLUST_kmer-", 
           kmer, "_", action, ".pdf"),
  height = 7, width = 12
)
par(mfrow = c(1,2))
mclust::plot.Mclust(custom_clust, what = "classification", main = TRUE)
plot(1:15, wss, 
  type="b", 
  xlab="Number of Clusters",
  ylab="Within groups sum of squares"
)
abline(v = 6, col = 'darkred', lwd = 2, lty = 2)
plots.saved <- dev.off()

p1 <- as_tibble(my_pca.clust.breaktype) %>% 
  ggplot(aes(x = PC1, y = PC2, col = as.factor(kmeans.clust))) + 
  geom_point() +
  ggrepel::geom_text_repel(
    aes(x = PC1, y = PC2, label = breakage_category),
    box.padding = 1.2, 
    max.overlaps = Inf,
    min.segment.length = unit(0.1, "lines"),
    segment.size = 0.5,
    segment.color = "grey",
    size = 2.5) + 
  theme_bw() + 
  theme_classic() + 
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  ) 

ggsave(
    filename = paste0("../figures/PCA/",
                      "Clustering_Breakage_kmer-", kmer, "_", 
                      action, ".pdf"), 
    plot = p1,
    width = 10, height = 8
)

# save clustering results
fwrite(
  my_pca.clust.breaktype, 
  paste0(
    "../data/kmertone_tracts/",
    "Clustering_PCA_kmer-", kmer, 
    "_", action, ".csv"
  )
)

################################################################################
if(kmer < 6){
  # Cluster PC1 and PC2 by kmer
  my_pca.clust.kmer <- data.table(
    kmer = attr(my_pca$x, "dimnames")[[1]],
    PC1 = my_pca$x[, 1],
    PC2 = my_pca$x[, 2]
  )

  # calculate within groups sum of squares
  wss <- (nrow(my_pca.clust.kmer)-1)*sum(apply(
    my_pca.clust.kmer[, -"kmer"],2,var
  ))

  set.seed(1234)
  for (i in 2:15){
    wss[i] <- sum(kmeans(
      my_pca.clust.kmer[, -"kmer"], centers = i)$withinss
    )
  }

  # finding optimal nr of clusters based on BIC
  set.seed(1234)
  d_clust <- Mclust(
    my_pca.clust.kmer[, -"kmer"], 
    G=1:15, 
    verbose = FALSE
  )
  m.best <- dim(d_clust$z)[2]

  # plot combined results
  pdf(paste0("../figures/PCA/",
            "Clustering_PCs_BreakageKmer_kmer-", 
            kmer, "_", action, ".pdf"),
    height = 8,
    width = 14
  )
  par(mfrow = c(1,2))
  plot.Mclust(d_clust, what = "classification", main = TRUE)
  plot(1:15, wss, 
    type="b", 
    xlab="Number of Clusters",
    ylab="Within groups sum of squares"
  )
  plots.saved <- dev.off()

  # label clusters
  my_pca.clust.kmer[, kmeans.clust := d_clust$classification]

  p1 <- as_tibble(my_pca.clust.kmer) %>% 
    ggplot(aes(x = PC1, y = PC2, col = as.factor(kmeans.clust))) + 
    geom_point() +
    ggrepel::geom_text_repel(
      aes(x = PC1, y = PC2, label = kmer),
      box.padding = 1.1, 
      max.overlaps = Inf,
      min.segment.length = unit(0.1, "lines"),
      segment.size = 0.5,
      segment.color = "grey",
      size = 2) + 
    theme_bw()

  ggsave(
      filename = paste0("../figures/PCA/",
                        "Clustering_BreakageKmer_kmer-", kmer, "_", 
                        action, ".pdf"), 
      plot = p1,
      width = 14,
      height = 12
  )
}