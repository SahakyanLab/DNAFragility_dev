suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(ggseqlogo)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(gridExtra)))
suppressPackageStartupMessages(suppressWarnings(library(cluster)))
suppressPackageStartupMessages(suppressWarnings(library(dendextend)))
suppressPackageStartupMessages(suppressWarnings(library(gplots)))
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

# data sourcing: size of samples
org_file <- fread("../../data/org_file.csv", showProgress = FALSE)
org_file[, exp := paste0(`Fragmentation type`, "/", `Experiment folder`)]
sample.sizes <- tibble(
    main = org_file$Category_general,
    source = org_file$Class) %>% 
    dplyr::group_by(main) %>% 
    dplyr::mutate(count = dplyr::n()) %>% 
    dplyr::mutate(main = ifelse(stringr::str_detect(
        string = main, 
        pattern = "Epigenome"), 
        "Epigenome_mark", main)) %>%
    dplyr::summarise(count = dplyr::n()) %>% 
    tibble::add_row(main = "QM_params", count = 3) %>% 
    dplyr::mutate(count = ifelse(main == "Epigenome_mark", count+1, count)) %>% 
    dplyr::filter(main != "G4_map") %>% 
    dplyr::arrange(count) %>% 
    dplyr::mutate(main = forcats::fct_inorder(main)) 

##########################################################################################
#' Correlation of z-scores between each dataset.
##########################################################################################
cor.mat <- lapply(c("zscore", "ratio"), function(action){
    res <- pbapply::pblapply(c(2,4,6,8), function(kmer){
        # import datasets
        files <- list.files(
            path = "../data/kmertone",
            pattern = paste0("score_", kmer, "-mers.csv"),
            full.names = TRUE,
            recursive = TRUE
        )
        files <- files[grepl(pattern = "00_Breakage", x = files)]

        # # only filter for TK6 cell-line
        # files <- files[stringr::str_detect(
        #     string = files, pattern = "TK6_Top2_mediated_DSBs/ETO_"
        # )]

        exp.labels <- stringr::str_remove_all(
            string = files,
            pattern = paste0("../data/kmertone/|/score_", 
                                kmer, "-mers.csv")
        )
        org.file <- fread("../../data/org_file.csv")
        org.file[, exp := paste0(`Fragmentation type`, "/", `Experiment folder`)]
        breaktype.labels <- org.file$Category_Main[match(exp.labels, org.file$exp)]

        # get correlation density plot of experiments from same type
        df <- tibble(
            file.name = files,
            exp = exp.labels,
            bp.type = breaktype.labels) %>% 
            dplyr::group_by(bp.type) %>% 
            dplyr::mutate(count = dplyr::n()) %>% 
            dplyr::filter(count > 1) %>% 
            dplyr::select(-count) %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate(
                exp = stringr::str_remove_all(
                    string = exp,
                    pattern = "00_Breakage/"
                )
            )
            
        files.df <- df %>% 
            dplyr::group_split(bp.type)

        cor.mat <- lapply(files.df, function(f){
            res <- lapply(f$file.name, function(x){
                temp <- switch(action,
                    "zscore" = {
                        fread(x, select = "z")
                    },
                    "ratio" = {
                        dt <- fread(x, select = c("case", "control"))
                        norm.case <- dt$case/sum(dt$case, na.rm = TRUE)
                        norm.control <- dt$control/sum(dt$control, na.rm = TRUE)
                        ratio <- norm.case/norm.control
                    }
                )
                return(temp)
            })
            res <- do.call(cbind, res)
            colnames(res) <- paste0("V", 1:ncol(res))

            # remove any values of inf and NAs
            rows.to.keep <- which(is.finite(rowSums(res)))
            res <- res[rows.to.keep, ]

            # calculate correlation coefficients
            cor.coef <- cor(
                res, 
                use = "everything", 
                method = "pearson"
            )
            cor.coef <- as.numeric(cor.coef[1,])

            # keep non 1 correlation coefficients
            cor.coef <- cor.coef[cor.coef != 1]    

            return(list(cor.coef, rep(unique(f$bp.type), length(cor.coef))))
        })

        cors <- unlist(sapply(cor.mat, `[[`, 1))
        exps <- unlist(sapply(cor.mat, `[[`, 2))
        res <- data.table(exp = exps, cors = cors, kmer = kmer, action = action)

        return(res)
    })
    res <- rbindlist(res)
    return(res)
})
cor.mat.original <- rbindlist(cor.mat)

cor.mat <- as_tibble(cor.mat.original) %>% 
    dplyr::mutate(
        action = ifelse(action == "ratio", 
            "Probability Ratio", "Z-score"),
        action = factor(action, levels = c("Z-score", "Probability Ratio")),
        kmer = as.factor(kmer))
    
p1 <- cor.mat %>% 
    ggplot(aes(x = cors, fill = kmer)) + 
    geom_density(alpha = 0.4) + 
    facet_wrap(vars(action), nrow = 1) + 
    coord_cartesian(xlim = c(0, 1)) + 
    scale_x_continuous(limits = c(0, NA)) +
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 20)) +
    labs(
        x = "Pearson correlation coefficient",
        y = "Density",
        title = paste0("Pearson correlation coefficients \n",
                       "between datasets of the same breakage type.")
    )

ggsave(
    filename = "../figures/correlation_maps/All_kmers.pdf",
    dpi = 600,
    plot = p1,
    height = 7,
    width = 12
)

##########################################################################################
#' Which enriched/depleted k-mers are highly pos/neg correlated across breakage types
##########################################################################################
# vars
args <- commandArgs(trailingOnly = TRUE)
kmer <- as.numeric(args[2])
action <- as.character(args[3])
only_breaks <- as.logical(args[4])

org_file <- fread("../../data/org_file.csv", showProgress = FALSE)
org_file[, exp := paste0(`Fragmentation type`, "/", `Experiment folder`)]

get_data <- function(kmer, action, only_breaks){
    # import processed kmertone data
    data.set <- fread(
        file = paste0(
            "../data/kmertone/QueryTable/QueryTable_kmer-",
            kmer, "_", action, ".csv"
        ),
        showProgress = FALSE
    )

    # subset only for DNA breakage datasets
    if(only_breaks){
        inds <- match(data.set$category, org_file$Category_Main)
        labels <- org_file$Class[inds]
        data.set <- data.set[which(labels == "DNA_Breakage"),]
    }

    # drop infinite values.
    # plus 1 because category is in the first column. 
    check.cols.with.inf <- which(is.infinite(colSums(data.set[,-"category"])))+1
    if(length(check.cols.with.inf) > 0) data.set <- data.set[, -check.cols.with.inf]

    return(data.set)
}

# combine results
df.zscore <- get_data(kmer=kmer, action=action, only_breaks=only_breaks)

# clustering requires any column to have 0 NAs. 
# Thus, remove rows of the column with least number of NAs.
check.nas <- apply(is.na(df.zscore[,-"category"]), 2, sum)
col.to.clean <- names(check.nas[which.min(check.nas)])
rows.to.keep <- which(!is.na(df.zscore[, ..col.to.clean]))
df.zscore <- df.zscore[rows.to.keep, ]

# normalise by row (experiment)
df.norm <- apply(df.zscore[, -"category"], 1, scale, 
                 center = TRUE, scale = TRUE)
df.norm <- t(df.norm)
rownames(df.norm) <- df.zscore$category
colnames(df.norm) <- colnames(df.zscore[, -"category"])
df.norm <- df.norm[,which(is.finite(colMeans(df.norm)))]

# cluster both by experiments and k-mers
# by: row (experiment)
df.hclust <- hclust(dist(df.norm), method = "ward.D2")
df.dendro <- as.dendrogram(df.hclust)

# by: column (k-mers)
df.kmers <- t(df.norm)
df.dendro.kmers <- as.dendrogram(hclust(
    dist(df.kmers), method = "ward.D2"
))

cols <- 10
color.scheme <- rev(brewer.pal(cols, "RdBu"))
CI <- mean(df.norm, na.rm = TRUE)+c(-1,1)*1.96*sd(df.norm, na.rm = TRUE)
# pdf("../figures/kmertone_tracts/for_paper.pdf", height = 10, width = 12)
png(paste0(
        "../figures/kmertone_tracts/for_paper_", 
        "kmer-", kmer, "_", action, ".png" 
    ),
    height = 10, width = 12, 
    units = "in", res = 300
)
par(oma=c(0,0,0,15))
out <- heatmap.2(
    df.norm,
    offsetRow = 0,
    offsetCol = 0,
    Rowv = df.dendro,
    Colv = df.dendro.kmers,
    dendrogram = "both",
    revC = FALSE,
    trace = "none",
    density.info = "none",
    margins = c(1,19),

    # labels
    key = FALSE,
    symkey = TRUE,
    keysize = 1.2,
    key.xlab = "Scaled z-scores",
    cexRow = 0.7,
    labRow = rownames(df.norm),
    labCol = "",               

    # colours
    col = color.scheme, 
    breaks = seq(CI[1], CI[2], length.out = cols+1),
    notecol = "black"
)
plot.save <- dev.off()

# colour dendrogram by cluster
ks <- 30
kmer.dend <- out$colDendrogram
coldend.clusters <- cutree(kmer.dend, k = ks)
coldend.clusters <- coldend.clusters[order.dendrogram(kmer.dend)]
coldend.clusters <- data.table(
    kmers = names(coldend.clusters),
    clusters = unname(coldend.clusters),
    hex.col = rainbow(n = ks)[unname(coldend.clusters)]
)
labels_colors(kmer.dend) <- coldend.clusters$hex.col 

# plot k-mer clusterings (columns)
png(
    paste0("../figures/kmertone_tracts/", kmer, 
           "-mer_", action, "-hclusteringONLY_byKmers.png"),
    height = 12, 
    width = 80,
    units = "in", res = 600
)
par(mar = c(5,3,5,0))
plot(kmer.dend %>% set("labels_cex", c(0.2, 0.2)))
plotted <- dev.off() 

################################################################################################################
# explore highly enriched/depleted k-mers from dendrogram/heatmap
# from top to bottom row, left to right.
cfDNA_breast_cancer <- df.norm[grepl("Breast", rownames(df.norm)),]
cfDNA_breast_cancer_kmers <- sort(cfDNA_breast_cancer, decreasing = TRUE)
most_fragile_kmers <- names(head(cfDNA_breast_cancer_kmers))
least_fragile_kmers <- names(tail(cfDNA_breast_cancer_kmers))

kmer_clusters <- c(
    coldend.clusters$clusters[
        coldend.clusters$kmers == least_fragile_kmers[1]
    ],
    coldend.clusters$clusters[
        coldend.clusters$kmers == most_fragile_kmers[1]
    ]
)

plots <- lapply(c(kmer_clusters), function(this.cluster){
    # extract subset of kmers from cluster
    kmer.subset <- coldend.clusters[clusters == this.cluster]$kmers

    # # get sequence logos
    # rev.comp <- Biostrings::reverseComplement(Biostrings::DNAStringSet(kmer.subset))
    # kmer.subset <- c(kmer.subset, paste(rev.comp))
    p1 <- suppressWarnings(ggseqlogo(kmer.subset)) + 
        coord_cartesian(ylim = c(0,2)) + 
        labs(y = "") + 
        theme(
            axis.line.y = element_line(colour = "black"),
            axis.text.x = element_blank(),
            # axis.text.y = element_text(size = 30) 
            axis.text.y = element_blank()
        )
    return(p1)
})

pdf(
    paste0("../figures/kmertone_tracts/", kmer, 
           "-mer_", action, "-hclusteringONLY_byKmers_SeqLogos.pdf"),
    height = 4, 
    width = 10
)
do.call(gridExtra::grid.arrange, c(plots, nrow = 1))
plotted <- dev.off()
################################################################################################################

# colour by breakage type
inds <- match(rownames(df.norm), org_file$Category_Main)
category.labels <- org_file$Category_Main[inds]
category.general.labels <- org_file$Category_general[inds]
category.general.hex <- org_file$Category_Hex[inds]
col.df <- data.table(
    exp = rownames(df.norm),
    category = category.labels,
    breakage.type = category.general.labels,
    hex.cols = category.general.hex
)

# evaluate optimal number of clusters
suppressPackageStartupMessages(library(factoextra))
max.ks <- 30
p1 <- fviz_nbclust(df.norm, hcut, method = "wss", 
                   k.max = max.ks, linecolor = "#6D7997") + 
    theme(text = element_text(size = 7)) + 
    coord_cartesian(xlim = c(1,max.ks))
p2 <- fviz_nbclust(df.norm, hcut, method = "silhouette", 
                   k.max = max.ks, linecolor = "#6D7997") + 
    coord_cartesian(xlim = c(1,max.ks)) + 
    theme(text = element_text(size = 7)) + 
    labs(title = "", x = "")

pdf(paste0(
    "../figures/kmertone_tracts/", kmer, 
    "-mer_", action, "-hclusteringONLY_byExp_eval.pdf"),
    height = 8, width = 12
)
gridExtra::grid.arrange(p1, p2, ncol = 2)
plot.saved <- dev.off()

# no strong convergence, hence will cluster based on the 
# 2 principal component results first
pca_clust <- fread(paste0(
    "../data/kmertone_tracts/",
    "Clustering_PCA_kmer-", kmer, 
    "_", action, ".csv"
))
max.ks <- max(pca_clust$kmeans.clust)
# match(rownames(df.norm), pca_clust$breakage_category)

tree.clusters <- cutree(df.dendro, k = max.ks)
col.df[, clusters := unname(tree.clusters)[match(
    col.df$exp, names(tree.clusters))]]

# order by dendogram
col.df <- col.df[order(match(
    exp, col.df$exp[order.dendrogram(df.dendro)]
))]
labels_colors(df.dendro) <- col.df$hex.cols

png(
    paste0("../figures/kmertone_tracts/", kmer, 
           "-mer_", action, "-hclusteringONLY_byExp_BlackLabels.png"),
    height = 10, 
    width = 9,
    units = "in", res = 600
)
par(mar = c(4,1,1,19))
plot(
    df.dendro %>% set("labels_cex", c(0.9, 0.9)),
    horiz = TRUE
)
plot.save <- dev.off()

png(
    paste0("../figures/kmertone_tracts/", kmer, 
           "-mer_", action, "-hclusteringONLY_byExp_withRECT.png"),
    height = 10, 
    width = 9,
    units = "in", res = 600
)
par(mar = c(4,1,1,19))
plot(
    df.dendro %>% set("labels_cex", c(0.9, 0.9)),
    horiz = TRUE
)
rect.dendrogram(df.dendro, k=max.ks, horiz=TRUE)
plot.save <- dev.off()

# 2. calculate correlation matrix between all k-mers and all clusters
# first, need to reorder numbers based on rank
ordered.factor <- factor(
    col.df$clusters, 
    levels = rev(unique(col.df$clusters)), 
    ordered = TRUE
)
decreasing.order <- as.numeric(ordered.factor)
col.df$clusters <- decreasing.order
setorder(col.df, clusters)

# average k-meric enrichment/depletion by cluster
res <- lapply(unique(col.df$clusters), function(x){
    df.norm.filtered.a <- df.norm[match(
        col.df$category[col.df$clusters == x], rownames(df.norm)
    ),]
    df.norm.filtered.a <- colMeans(df.norm.filtered.a, na.rm = TRUE)
})

# comparing cluster 1 vs. cluster 2 is the same as cluster 2 vs. cluster 1.
# Hence, only keep unique comparisons. 
cor.mat.ind <- expand.grid(
    x = unique(col.df$clusters),
    y = unique(col.df$clusters)
)
cor.mat.ind <- cor.mat.ind[cor.mat.ind[,1] != cor.mat.ind[,2],]
cor.mat.ind <- as_tibble(cor.mat.ind) %>% dplyr::group_split(y)
cor.mat.ind <- head(cor.mat.ind, n = -1)
cor.mat.ind <- lapply(1:length(cor.mat.ind), function(x){
    as.data.frame(cor.mat.ind[[x]][x:nrow(cor.mat.ind[[x]]),])
})
cor.mat.ind <- do.call(rbind, cor.mat.ind)
rownames(cor.mat.ind) <- NULL

# get hybridisation energies
if(kmer >= 6){
    all.kmers <- do.call(
        data.table::CJ, 
        rep(list(c("A", "C", "G", "T")), kmer)
    )
    kmer_list <- all.kmers[, do.call(paste0, .SD)]
    kmer_ref <- data.table(
        'kmer' = kmer_list,
        'rev.comp' = as.character(
            Biostrings::reverseComplement(
                Biostrings::DNAStringSet(kmer_list)
            )
        )
    )
    kmer_ref[, `:=`(cond = ifelse(seq(1:nrow(.SD)) < match(kmer, rev.comp), 
    TRUE, ifelse(kmer == rev.comp, TRUE, FALSE)))]
    kmer_ref <- kmer_ref[cond == TRUE, .(kmer, rev.comp)]

    # get delta heat of formation
    get_qm_params = function(kmer){
        # only keep heat of formation and ionisation potential differences
        df <- fread("../../data/04_QM_parameters/6mer/Raw/denergy.txt")
        cols.to.keep <- c(1, which(grepl(pattern = "B_ds.dEhof", x = colnames(df))))
        df <- df[, ..cols.to.keep]

        if(kmer == 8){
            # if breakage k-mer > QM kmers, then use a sliding window approach
            query.kmers <- kmer_list
            len <- 3

            # rolling window substring
            list.kmer.vals <- lapply(1:len, function(x){
                first.kmer <- substring(
                    text = query.kmers,
                    first = x, 
                    last = x+5
                )
                first.kmer.ind <- match(first.kmer, df$seq)
                first.kmer.vals <- df[first.kmer.ind, -"seq"]
            })

            # average all
            list.kmer.vals <- list.kmer.vals[[1]]+list.kmer.vals[[2]]+list.kmer.vals[[3]]
            df <- list.kmer.vals/len
            df[, seq := query.kmers]
        }
        return(df)
    }

    source("./01_HybridisationGfree/GEN_getKmerHybridisationGs.R")
    df.hybridenergies <- getKmerHybridisationGs(k = kmer, plot = FALSE)
    
    df.dehof <- get_qm_params(kmer = kmer)
    df.dehof.kmers <- df.dehof$seq

    # normalise to similar scales as GFrees
    df.dehof <- data.table(
        kmer = df.dehof.kmers, 
        dehof = scale(df.dehof$B_ds.dEhof_ds_ss)[,1]
    )

    # first first lexicologically occurring k-mer
    rows_to_keep <- match(kmer_ref$kmer, df.dehof$kmer)
    df.dehof <- df.dehof[rows_to_keep,]
    df.hybridenergies <- df.hybridenergies[rows_to_keep,]

    # extract the most different k-mers
    df.res <- lapply(1:nrow(cor.mat.ind), function(x){
        ind.x <- cor.mat.ind[x,"x"]
        ind.y <- cor.mat.ind[x,"y"]
        res.both <- as.data.frame(cbind(x = res[[ind.x]], y = res[[ind.y]]))
        res.both$xy.diff <- res.both$y-res.both$x
        res.both <- res.both[order(res.both$xy.diff, decreasing = TRUE),]

        # extract top 1% extreme k-mers, i.e. k-mers that are most different
        top.05.perc.each.side <- ceiling(nrow(res.both)*0.005)
        top.1.perc <- rbind(
            res.both[1:top.05.perc.each.side,],
            res.both[(nrow(res.both)-top.05.perc.each.side+1):nrow(res.both),]
        )

        top.1.perc <- top.1.perc %>% 
            tibble::rownames_to_column("kmer") %>% 
            as_tibble() %>% 
            dplyr::select(kmer, xy.diff) %>% 
            dplyr::arrange(xy.diff) %>% 
            dplyr::mutate(kmer = forcats::fct_inorder(kmer)) %>% 
            dplyr::mutate(
                cluster = paste0("Cluster: ", ind.y, " & ", ind.x),
            )
        top.1.perc <- as.data.table(top.1.perc)

        # plot boxplots against hybridisation energies
        top.1.perc[, `:=`(
            hybridisation.Gfrees = df.hybridenergies$Gfrees[match(
                top.1.perc$kmer, df.hybridenergies$kmers
            )],
            dEhof = df.dehof$dehof[match(top.1.perc$kmer, df.dehof$kmer)],
            enriched = ifelse(xy.diff >= 0, "Enriched", "Depleted"),
            col.hex = ifelse(xy.diff >= 0, "#b2182b", "#2166ac")
        )]

        # plot extreme values
        res.both <- res.both %>% 
            tibble::rownames_to_column("kmers") %>% 
            as_tibble() %>% 
            dplyr::mutate(
                extreme.kmer = ifelse(
                    row_number() <= top.05.perc.each.side | 
                    row_number() >= (dplyr::n()-top.05.perc.each.side+1), 
                    TRUE, FALSE),
                col.hex = ifelse(
                    row_number() <= top.05.perc.each.side, "#b2182b", 
                    ifelse(row_number() >= (dplyr::n()-top.05.perc.each.side+1), "#2166ac",
                    "#dddddd"))
            )
        res.both.labels <- res.both %>% dplyr::filter(extreme.kmer == TRUE) 

        p1 <- res.both %>% 
            ggplot(aes(x = xy.diff, fill = col.hex)) + 
            geom_histogram(bins = 100) + 
            # ggrepel::geom_text_repel(
            #     data = res.both.labels, 
            #     aes(x = xy.diff, y = 0, label = kmers),
            #     box.padding = 1.5, 
            #     max.overlaps = Inf,
            #     min.segment.length = unit(0.01, "lines"),
            #     segment.color = "grey",
            #     segment.size = 0.1,
            #     size = 4
            # ) + 
            theme_bw() +
            theme(
                panel.border = element_blank(),
                panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.position = "none",
                text = element_text(size = 30),
                plot.title = element_text(size = 20)) + 
            scale_fill_identity() + 
            labs(
                title = paste0("Cluster: ", ind.y, " - ", ind.x),
                x = "", 
                y = ""
            )

        # example plot for publication figure: histogram
        if(all(c(1,2) %in% c(ind.x, ind.y))){
            p1 <- p1 + 
                coord_cartesian(xlim = c(-5,5))
            ggsave(
                filename = paste0(
                    "../figures/kmertone_tracts/", 
                    kmer, "-mer_", action, 
                    "-mostextremekmers_HISTOGRAM_HybridFE_example_1_2.pdf"
                ),
                plot = p1,
                height = 5, width = 5
            )
        }

        return(list(top.1.perc, p1))
    })
    df.most.diff.kmers <- rbindlist(sapply(df.res, `[`, 1))

    df_most_diff_kmers <- as_tibble(df.most.diff.kmers) %>% 
        dplyr::mutate(cluster = stringr::str_replace(
            string = cluster, pattern = "&", replacement = "-"
        )) %>%
        dplyr::mutate(cluster = stringr::str_remove(
            string = cluster, pattern = "Cluster: "
        )) %>% 
        dplyr::mutate(cluster = forcats::fct_inorder(cluster))

    # as_tibble(df.res[[1]][[1]])
    #     dplyr::mutate(cluster = stringr::str_replace(
    #         string = cluster, pattern = "&", replacement = "-"
    #     )) %>%
    #     dplyr::mutate(cluster = stringr::str_remove(
    #         string = cluster, pattern = "Cluster: "
    #     )) %>% 
    #     dplyr::mutate(cluster = forcats::fct_inorder(cluster))

    plot_most_diff_kmers <- function(
        dat, col_of_int, filter_for_cluster, 
        x_factor, y_factor, scatter_size = 1, plot_count = FALSE){
        var_of_int <- dat %>% 
            dplyr::select({{col_of_int}}) %>% 
            # dplyr::select("hybridisation.Gfrees") %>%
            dplyr::pull(.)

        if(!is.null(filter_for_cluster)){
            dat <- dat %>% 
                dplyr::filter(cluster == filter_for_cluster)
        }

        y_range <- x_range <- NA
        if(!is.null(x_factor)) x_range <- min(var_of_int)*x_factor
        if(!is.null(y_factor)) y_range <- max(var_of_int)*y_factor

        if(plot_count){
            p1 <- dat %>%
                # ggplot(aes(x = {{col_of_int}})) + 
                ggplot(aes(x = hybridisation.Gfrees)) + 
                geom_histogram(col = "grey80") + 
                facet_wrap(vars(cluster), nrow = 1) +
                coord_cartesian(ylim = c(x_range, y_range)) +
                theme_bw() + 
                theme(
                    text = element_text(size = 20),
                    axis.text.x = element_text(size = 20),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    strip.background = element_blank(),
                    legend.position = "none",
                    strip.text = element_text(size = 15)
                ) +
                scale_colour_identity() + 
                scale_x_discrete(labels = c(
                    "Depleted" = "-",
                    "Enriched" = "+")
                ) + 
                labs(
                    x = "",
                    y = ""
                )

        } else {
            p1 <- dat %>%
                ggplot(aes(x = enriched, y = {{col_of_int}})) + 
                geom_boxplot(outlier.shape = NA) + 
                geom_jitter(
                    aes(col = col.hex), size = scatter_size,
                    position = position_jitter(width = 0.1, height = 0.1)) +
                ggsignif::geom_signif(
                    comparisons = list(c("Enriched", "Depleted")),
                    map_signif_level = TRUE,
                    test = "t.test",
                    textsize = 5,
                    vjust = -0.1,
                    margin_top = 0.1
                ) +
                facet_wrap(vars(cluster), nrow = 1) +
                coord_cartesian(ylim = c(x_range, y_range)) +
                theme_bw() + 
                theme(
                    text = element_text(size = 20),
                    axis.text.x = element_text(size = 20),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    strip.background = element_blank(),
                    legend.position = "none",
                    strip.text = element_text(size = 15)
                ) +
                scale_colour_identity() + 
                scale_x_discrete(labels = c(
                    "Depleted" = "-",
                    "Enriched" = "+")
                ) + 
                labs(
                    x = "",
                    y = ""
                )
        }

        return(p1)
    }

    dir.create(
        path = "../data/kmertone_tracts/",
        showWarnings = FALSE
    )
    fwrite(
        df.most.diff.kmers,
        file = paste0("../data/kmertone_tracts/", kmer, 
                    "-mer_most_different_kmers.csv")
    )

    # plot: hybridisation energy and delta heat of formation
    suppressPackageStartupMessages(library(ggsignif))
    ggsave(
        filename = paste0(
            "../figures/kmertone_tracts/", 
            kmer, "-mer_", action, 
            "-mostextremekmers_BOXPLOTS_HybridFE.pdf"
        ),
        plot = plot_most_diff_kmers(
            dat = df_most_diff_kmers,
            col_of_int = hybridisation.Gfrees,
            filter_for_cluster = NULL,
            x_factor = 1.1,
            y_factor = 0.3
        ),
        height = 4, width = 16
        # dpi = 600
    )

    # plot: delta heat of formation
    ggsave(
        filename = paste0(
            "../figures/kmertone_tracts/", 
            kmer, "-mer_", action, 
            "-mostextremekmers_BOXPLOTS_dEhof.pdf"
        ),
        plot = plot_most_diff_kmers(
            dat = df_most_diff_kmers,
            col_of_int = dEhof,
            filter_for_cluster = NULL,
            x_factor = 1.1,
            y_factor = 1.45
        ),
        height = 4, width = 16
        # dpi = 600
    )

    # example plot for publication figure: boxplot
    ggsave(
        filename = paste0(
            "../figures/kmertone_tracts/", 
            kmer, "-mer_", action, 
            "-mostextremekmers_BOXPLOTS_HybridFE_example_1_2.pdf"
        ),
        plot = plot_most_diff_kmers(
            dat = df_most_diff_kmers,
            col_of_int = hybridisation.Gfrees,
            filter_for_cluster = "1 - 2",
            x_factor = NULL,
            y_factor = 0.5,
            scatter_size = 3
        ),
        height = 5, width = 4
        # dpi = 600
    )

    # extract p-values
    example_test <- df_most_diff_kmers %>% 
        dplyr::filter(cluster == "1 - 2") %>% 
        dplyr::select(enriched, hybridisation.Gfrees) %>% 
        dplyr::group_split(enriched) %>% 
        purrr::set_names(purrr::map(., ~ unique(.x$enriched)[1]))
    example_test <- t.test(
        x = example_test$Enriched$hybridisation.Gfrees, 
        y = example_test$Depleted$hybridisation.Gfrees,
        alternative = "two.sided"
    )
    # example_test$p.value = 5.912104e-37

    # example plots for publication figure
    ggsave(
        filename = paste0(
            "../figures/kmertone_tracts/", 
            kmer, "-mer_", action, 
            "-mostextremekmers_BOXPLOTS_HybridFE_COUNT_example_1_2.pdf"
        ),
        plot = plot_most_diff_kmers(
            dat = df_most_diff_kmers,
            col_of_int = hybridisation.Gfrees,
            filter_for_cluster = "1 - 2",
            x_factor = NULL,
            y_factor = 0.5,
            scatter_size = 3,
            plot_counts = TRUE
        ),
        height = 5, width = 4
        # dpi = 600
    )

    # # plot: histogram with annotated extreme k-mer
    # png(
    #     paste0("../figures/kmertone_tracts/", kmer, 
    #            "-mer_", action, "-mostextremekmers.png"),
    #     height = 25, 
    #     width = 60,
    #     units = "in", res = 600
    # )
    # do.call(gridExtra::grid.arrange, c(sapply(df.res, `[`, 2), ncol = 10))
    # plotted <- dev.off()
}

# ##########################################################################################
# #' K-meric z-scores vs. probabilty ratios across breakage types
# ##########################################################################################
# # combine results
# df.zscore <- get_data(kmer=kmer, action="zscore", only_breaks=TRUE)
# df.ratio <- get_data(kmer=kmer, action="ratio", only_breaks=TRUE)

# # combine by kmers and category
# df.zscore <- df.zscore %>% tidyr::pivot_longer(
#     cols = !category,
#     names_to = 'kmer',
#     values_to = 'value'
# )
# df.ratio <- df.ratio %>% tidyr::pivot_longer(
#     cols = !category,
#     names_to = 'kmer',
#     values_to = 'value'
# )

# df <- inner_join(
#     x = df.zscore, 
#     y = df.ratio,
#     by = c("category", "kmer")) %>% 
#     dplyr::rename_with(~c(
#         "category", "kmer", 
#         "zscore", "ratio"
#     ))

# # volcano plots
# df <- df %>% 
#     dplyr::mutate(
#         zscore = abs(zscore), ratio = log2(ratio),
#         enriched = ifelse(ratio >= 0, TRUE, FALSE),
#         label.hex = ifelse(enriched, "#bb357e", "#4592c8")
#     )

# # drop infinite values
# check.rows.with.inf <- which(is.infinite(rowSums(df["ratio"])))
# if(length(check.rows.with.inf) > 0) df <- df[-check.rows.with.inf,]

# # drop nas
# df <- df %>% tidyr::drop_na()

# # assign overarching breakage category naming
# inds <- match(df$category, org_file$Category_Main)
# labels <- org_file$Category_general[inds]
# df <- df %>% dplyr::mutate(b.type = labels)

# # plot each breakage type separately
# df.split <- df %>% 
#     dplyr::group_split(b.type)

# df.bio.split <- df.split[[1]] %>% 
#     dplyr::mutate(is.wrn = ifelse(stringr::str_detect(
#         string = category, pattern = "WRN_loss"), 
#         TRUE, FALSE)) %>% 
#     dplyr::group_split(is.wrn)

# df.split <- df.split[-1]
# df.split <- c(df.split, df.bio.split)

# if(kmer == 8){
#     for(i in 1:length(df.split)){
#         df.split[[i]]$category <- stringr::str_replace_all(
#             string = df.split[[i]]$category,
#             pattern = "_",
#             replacement = " "
#         )

#         p1 <- df.split[[i]] %>% 
#             ggplot(aes(x = ratio, y = zscore)) +
#             geom_point(
#                 col = df.split[[i]]$label.hex, 
#                 alpha = 0.2
#             ) +    
#             facet_wrap(
#                 vars(category), 
#                 labeller = label_wrap_gen(), 
#                 scales = "free_y",
#                 ncol = 7
#             ) +
#             theme_bw() + 
#             theme_classic() + 
#             theme(text = element_text(size = 13)) + 
#             labs(
#                 x = "Log2 probability ratio",
#                 y = "Absolute z-scores"
#             )

#         dir.create(
#             path = "../figures/volcano_plots",
#             showWarnings = FALSE,
#             recursive = TRUE
#         )

#         filename <- paste0(
#             "../figures/volcano_plots/", 
#             unique(df.split[[i]]$b.type),
#             "_kmer-", kmer, ".png"
#         )

#         switch(unique(df.split[[i]]$b.type),
#             "Mechanical" = ggsave(filename = filename, plot = p1, height = 4, width = 7),
#             "Natural Fragmentation" = ggsave(filename = filename, plot = p1, height = 4, width = 12),
#             "Cell free DNA" = ggsave(filename = filename, plot = p1, height = 7, width = 14), 
#             "Enzymatic" = ggsave(filename = filename, plot = p1, height = 4, width = 12),
#             "Biological" = ggsave(filename = filename, plot = p1, height = 8, width = 16)
#         )
#     }
# }

##########################################################################################
#' Plot sequence logos for each experiment
##########################################################################################
# df.unique.cat <- df %>% 
#     dplyr::distinct(category, .keep_all = TRUE) %>% 
#     dplyr::select(category, b.type)

# df %>% 
#     dplyr::filter(b.type == "Biological") %>% 
#     dplyr::filter(category == "neuronal_cells") %>% 
#     dplyr::arrange(desc(zscore))

# # seq logos of top 5% enriched and top 5% depleted kmers
# for(i in 1:nrow(df.unique.cat)){
#     # filter data set
#     df.filtered <- dplyr::filter(df, category == df.unique.cat$category[i])

#     # top 5% of all
#     proportion <- 0.05
#     to.subsample <- ceiling(nrow(df.filtered)*proportion)
#     if(nrow(df.filtered)*0.05 < 5) to.subsample <- nrow(df.filtered)

#     # top 5% most enriched k-mers
#     df.enriched <- df.filtered %>% 
#         dplyr::arrange(desc(zscore), desc(ratio))
#     enriched.kmers <- df.enriched %>% 
#         dplyr::slice(1:to.subsample, with_ties = TRUE) %>% 
#         dplyr::pull(kmer)
#     enriched.kmers <- c(enriched.kmers, as.character(
#         Biostrings::reverseComplement(Biostrings::DNAStringSet(enriched.kmers)
#     )))

#     # plot seqlogos for enriched kmers
#     p1 <- ggseqlogo(enriched.kmers, method = "bits") + 
#         ylim(0, 2) +
#         labs(
#             title = "Top 5% enriched kmers",
#             subtitle = df.unique.cat$category[i]) %>% 
#         suppressWarnings()

#     # get top 5% depleted by probability ratios kmers
#     df.depleted <- df.filtered %>% 
#         dplyr::arrange(desc(zscore), ratio)
#     depleted.kmers <- df.depleted %>% 
#         dplyr::slice(1:to.subsample, with_ties = TRUE) %>% 
#         dplyr::pull(kmer)

#     # plot seqlogos for depleted kmers
#     p2 <- ggseqlogo(depleted.kmers, method = "bits") + 
#         ylim(0, 2) +
#         labs(title = "Top 5% depleted kmers") %>% 
#         suppressWarnings()

#     # save plots
#     dir.create(
#         path = paste0("../figures/seqlogos/kmer-", kmer, "/",
#                       df.unique.cat$b.type[i]),
#         showWarnings = FALSE,
#         recursive = TRUE
#     )
#     pdf(paste0("../figures/seqlogos/", 
#                "/kmer-", kmer, "/", 
#                df.unique.cat$b.type[i], "/",
#                df.unique.cat$category[i], ".pdf"),
#         height = 7,
#         width = 12
#     )
#     grid.arrange(p1, p2, ncol = 2)
#     plots.saved <- dev.off()
# }

##########################################################################################
#' Correlation heatmap between k-meric enrichment/depletion of each breakage types
##########################################################################################
suppressPackageStartupMessages(suppressWarnings(library(dendextend)))
org_file <- fread("../../data/org_file.csv", showProgress = FALSE)
org_file[, exp := paste0(`Fragmentation type`, "/", `Experiment folder`)]

df.zscore <- get_data(kmer=kmer, action="zscore", only_breaks=FALSE)
df.ratio <- get_data(kmer=kmer, action="ratio", only_breaks=FALSE)

for(action in c("zscore", "ratio")){
    df.metric <- switch(action,
        "zscore" = df.zscore,
        "ratio" = df.ratio,
    )

    # drop infinite values
    check.cols.with.inf <- which(is.infinite(colSums(df.metric[, -"category"])))+1
    if(length(check.cols.with.inf) > 0) df.metric <- df.metric[, -check.cols.with.inf]

    # drop nas
    df.metric <- df.metric %>% tidyr::drop_na()

    # standardise for correlation matrix
    df.metric.cat <- df.metric$category
    df.norm <- apply(df.metric[, -"category"], 2, scale, 
                    center = TRUE, scale = TRUE)
    rownames(df.norm) <- df.metric$category
    df.norm <- t(df.norm)

    # correlation heatmap with hierarchical clustering
    df.cor <- cor(
        df.norm, 
        use = "everything", 
        method = "pearson"
    )
    df.dist <- dist(df.cor)
    df.tree <- hclust(df.dist, method = "ward.D2")
    df.dend <- as.dendrogram(df.tree)
    color.scheme <- rev(brewer.pal(10, "RdBu"))

    # # reorder to match dendrogram's ordering
    # col.df <- data.table(exp = rownames(df.cor)[order.dendrogram(df.dend)])
    
    # colour by category
    inds <- match(rownames(df.cor), org_file$Category_Main)
    labels <- org_file$Category_general[inds]
    labels[which(is.na(labels))] <- "QM_parameters"

    label.hex <- org_file$Category_Hex[inds]
    label.hex[which(is.na(label.hex))] <- "#8F9CA7"

    col.df <- data.table(
        exp = rownames(df.cor),
        cols = label.hex
    )
    # col.df[, cols := label.hex]

    # # reorder to match dendrogram's ordering
    # rev.col.df <- copy(col.df)
    # rev.col.df[, ID := 1:nrow(rev.col.df)]
    # setorder(rev.col.df, -ID)
    # rev.col.df[, ID := NULL]

    if(only_breaks){
        dir.create(
            path = "../figures/correlation_maps",
            showWarnings = FALSE,
            recursive = TRUE
        )
        pdf(paste0("../figures/correlation_maps/", 
                    "/Hierarchical_clustering_of_correlations", 
                    "_byBreakageType_kmer-", kmer, 
                    "_", action, ".pdf"),
            height = 18,
            width = 20
        )
        heatmap.2(
            df.cor, 
            Rowv = df.dend,
            Colv = df.dend,
            dendrogram = "both", 
            revC = TRUE, 
            trace = "none", 
            density.info = "none",
            col = color.scheme, 
            key = TRUE, 
            key.xlab = "Correlation Coefficient",
            lhei = c(1,6), 
            lwid = c(1,6),
            margins = c(22,22),
            symkey = TRUE,
            cexRow = 1,
            cexCol = 1
            # colRow = col.df$cols,
            # labRow = col.df$exp,
            # colCol = col.df$cols
            # labCol = col.df$exp
        )
        plot.save <- dev.off()
    } else {
        pdf(paste0("../figures/correlation_maps/", 
                    "/Hierarchical_clustering_of_correlations", 
                    "_byAll_kmer-", kmer, 
                    "_", action, ".pdf"),
            height = 14,
            width = 16
        )
        heatmap.2(
            df.cor, 
            Rowv = df.dend,
            Colv = df.dend,
            dendrogram = "both", 
            revC = TRUE, 
            trace = "none", 
            density.info = "none",
            col = color.scheme, 
            key = TRUE, 
            key.xlab = "Correlation Coefficient",
            lhei = c(1,6), 
            lwid = c(1,6),
            margins = c(1,1),
            symkey = TRUE,
            cexRow = 1,
            cexCol = 1,
            labRow = "",
            labCol = ""
        )
        plot.save <- dev.off()

        pdf(paste0("../figures/correlation_maps/", 
                    "/Hierarchical_clustering_of_correlations", 
                    "_byAllNamed_kmer-", kmer, 
                    "_", action, ".pdf"),
            height = 42,
            width = 45
        )
        heatmap.2(
            df.cor, 
            Rowv = df.dend,
            Colv = df.dend,
            dendrogram = "both", 
            revC = TRUE, 
            trace = "none", 
            density.info = "none",
            col = color.scheme, 
            key = TRUE, 
            key.xlab = "Correlation Coefficient",
            lhei = c(1,6), 
            lwid = c(1,6),
            margins = c(22,22),
            symkey = TRUE,
            cexRow = 0.4,
            cexCol = 0.4
        )
        plot.save <- dev.off()
    }    
}

##########################################################################################
#' Create network graph in cytoscape from correlation matrix between all averaged datasets
##########################################################################################
# cytoscape 
suppressPackageStartupMessages(library(RCy3))
cytoscapeVersionInfo()

# list of app to install
cyto_app_toinstall <- c("clustermaker2", "enrichmentmap", "autoannotate")
cytoscape_version <- unlist(strsplit(cytoscapeVersionInfo()['cytoscapeVersion'],split = "\\."))
for(i in 1:length(cyto_app_toinstall)){
    print(i)
    commandsGET(paste("apps install app=\"", cyto_app_toinstall[i],"\"", sep=""))
}

only_breaks = FALSE
cor_threshold = 0.7
df.metric <- get_data(kmer=kmer, action="zscore", only_breaks=only_breaks)

# drop infinite values
check.cols.with.inf <- which(is.infinite(colSums(df.metric[, -"category"])))+1
if(length(check.cols.with.inf) > 0) df.metric <- df.metric[, -check.cols.with.inf]

# drop nas
df.metric <- df.metric %>% tidyr::drop_na()

# only keep DNA breakage and TFBS classes
inds <- match(df.metric$category, org_file$Category_Main)
labels <- org_file$Category_general[inds]
labels[which(is.na(labels))] <- "QM_parameters"
gen.class <- org_file$Class[inds]
gen.class[which(is.na(gen.class))] <- "QM_parameters"

# normalise by column because I want to keep the k-mer distribution
# per dataset in-tact
df.metric.cat <- df.metric$category
df.norm <- apply(df.metric[, -"category"], 2, scale, 
                 center = TRUE, scale = TRUE)
rownames(df.norm) <- df.metric$category

# correlation heatmap with hierarchical clustering
df.cor <- cor(
    t(df.norm), 
    use = "everything", 
    method = "pearson"
)

dir.create(path = "../data/correlations/", showWarnings = FALSE)
if(!only_breaks){
    # set the diagonal of matrix to zero - eliminate self-correlation
    df.cor[row(df.cor) == col(df.cor)] <- 0

    # set all correlations that are less than 0.9 to zero
    df.cor[which(df.cor <= cor_threshold & df.cor >= -cor_threshold)] <- 0

    # get rid of rows and columns that have no correlations with the above thresholds
    df.cor <- df.cor[which(rowSums(df.cor) != 0), which(colSums(df.cor) != 0)]

    correlation_filename <- paste0(getwd(), "/../data/correlations/",
                                   "correlation_map_all_", 
                                   kmer, "mer.txt")
} else {
    correlation_filename <- paste0(getwd(), "/../data/correlations/",
                                   "correlation_map_onlybreaks_", 
                                   kmer, "mer.txt")
}

# write out the correlation file
write.table(
    df.cor,  
    file = correlation_filename, 
    col.names = TRUE, 
    row.names = FALSE, 
    sep = "\t", 
    quote = FALSE
)

# generate network graph in cytoscape
amat_url <- "aMatReader/v1/import"
amat_params = list(
    files = list(correlation_filename),
    delimiter = "TAB",
    undirected = TRUE,
    ignoreZeros = TRUE,
    interactionName = "correlated with",
    rowNames = FALSE
)

response <- cyrestPOST(
    operation = amat_url, 
    body = amat_params, 
    base.url = "http://localhost:1234"
)
current_network_id <- response$data["suid"]

# colour by category
inds <- match(rownames(df.cor), org_file$Category_Main)
labels <- org_file$Category_general[inds]
labels[which(is.na(labels))] <- "QM_parameters"
label.hex <- org_file$Category_Hex[inds]
label.hex[which(is.na(label.hex))] <- "#8F9CA7"

node_class <- trimws(labels)
node_class_df <- data.frame(
    name = rownames(df.cor),
    node_class
)

# map new node attribute and categories to the network
loadTableData(
    node_class_df,
    table.key.column = "name",
    data.key.column = "name"  #default data.frame key is row.names
)

layoutNetwork(paste(
    "force-directed",
    "defaultSpringCoefficient=0.00003",
    "defaultSpringLength=70"
))

setNodeBorderColorDefault("#000000", style.name = "default")
setNodeLabelColorDefault("#000000", style.name = "default")
setNodeBorderWidthDefault(3, style.name = "default")
# setNodeShapeDefault('ellipse', style.name = "default")
setBackgroundColorDefault('#FFFFFF', style.name = "default")

# create colour map
col.df <- data.table(
    exp = rownames(df.cor),
    labels = labels,
    node.cols = label.hex,
    node.opacity = 255
)

# temporarily replace tfbs with higher opacities
node_opacity <- 255
col.df[, node.opacity := ifelse(node.cols == "#619B61", node_opacity, node.opacity)] # tfbs
col.df[, node.opacity := ifelse(node.cols == "#612940", node_opacity, node.opacity)] # epigenome faire-seq
col.df[, node.opacity := ifelse(node.cols == "#D15137", node_opacity, node.opacity)] # epigenome atac-seq
col.df[, node.opacity := ifelse(node.cols == "#E4D00A", node_opacity, node.opacity)] # epigenome chip-seq
col.df[, node.opacity := ifelse(node.cols == "#89CFF0", node_opacity, node.opacity)] # epigenome dnase-seq

# col.df[, node.opacity := 100]
setNodeFillOpacityMapping(
    table.column = "node_class",
    table.column.values = col.df$labels,
    opacities = col.df$node.opacity,
    mapping.type = "d",
    style.name = "default"   
)
setNodeBorderOpacityMapping(
    table.column = "node_class",
    table.column.values = col.df$labels,
    opacities = col.df$node.opacity,
    mapping.type = "d",
    style.name = "default" 
)
setNodeColorMapping(
    table.column = "node_class",
    table.column.values = col.df$labels,
    colors = col.df$node.cols,
    mapping.type = "d",
    style.name = "default"
)

# adjust label sizes if tfbs or epigenome
# col.df[, label.sizes := ifelse(grepl(pattern = "Epigenome|TFBS", x = labels), 0, 50)]
col.df[, label.sizes := 0] # uncomment if no labels to be plotted
setNodeFontSizeMapping(
    table.column = "node_class",
    table.column.values = col.df$labels,
    sizes = col.df$label.sizes,
    mapping.type = "d",
    style.name = "default" 
)

# colour edge based on positive or negative correlations
edge.table <- getTableColumns(
    table = "edge", 
    network = as.numeric(current_network_id)
)
col.with.corrs <- which(grepl(
    pattern = "correlation_map", 
    x = colnames(edge.table)
))
edge.table["edge_cols"] <- ifelse(edge.table[,col.with.corrs] >= 0, "#eb6979", "#93beea")
setEdgeColorMapping(
    table.column = basename(correlation_filename),
    table.column.values = edge.table[,col.with.corrs],
    colors = edge.table$edge_cols,
    mapping.type = "d",
    network = as.numeric(current_network_id),
    style.name = "default"
)

# node shape: circle
# node size: 200
# node transparency: 255
# edge transparency: 60
# edge width: 40

###########################################################################
# Highest absolute correlations of any TFBS with any type of breakage
###########################################################################
#' @description 
#' Passes the columns of interest into function and computes the 
#' correlation coefficient matrix.
#' @param dat filtered data table
#' @return filtered tibble data frame.
correlate.and.filter <- function(dat, rows.filter = "TFBS"){
    # correlation heatmap with hierarchical clustering
    df.cor <- cor(
        t(dat), 
        use = "everything", 
        method = "pearson"
    )

    df.cor <- as_tibble(df.cor) %>% 
        dplyr::mutate(
            exp = rownames(df.cor),
            .before = 1
        )

    # only keep breakage types in rows
    df.cor <- df.cor[which(!grepl(
        pattern = rows.filter, x = df.cor$exp)
    ),]

    # only keep tfbs in columns
    df.cor <- df.cor[c(1, which(grepl(
        pattern = rows.filter, x = colnames(df.cor)))
    )]
    df.cor <- df.cor %>% 
        tibble::column_to_rownames("exp") %>% 
        t() %>% 
        as.data.frame()

    # only keep highest absolute correlations
    df.cor[(df.cor < 0.7 & df.cor > -0.7)] <- NA

    # filter for finite values 
    df.cor.filtered <- df.cor %>% 
        tibble::rownames_to_column("tfbs") %>%
        as_tibble() %>%
        tidyr::pivot_longer(
            cols = !tfbs, 
            names_to = "DNA_breakage", 
            values_to = "corr") %>%
        dplyr::filter(is.finite(corr)) %>%
        tidyr::pivot_wider(
            names_from = "DNA_breakage", 
            values_from = "corr"
        )
    
    return(df.cor.filtered)
}

df.metric <- get_data(kmer=kmer, action="zscore", only_breaks=FALSE)

# drop infinite values
check.cols.with.inf <- which(is.infinite(colSums(df.metric[, -"category"])))+1
if(length(check.cols.with.inf) > 0) df.metric <- df.metric[, -check.cols.with.inf]

# drop nas
df.metric <- df.metric %>% tidyr::drop_na()

# only keep DNA breakage and TFBS classes
inds <- match(df.metric$category, org_file$Category_Main)
labels <- org_file$Category_general[inds]
labels[which(is.na(labels))] <- "QM_parameters"
gen.class <- org_file$Class[inds]
gen.class[which(is.na(gen.class))] <- "QM_parameters"

# normalise by column because I want to keep the k-mer distribution
# per dataset in-tact
df.metric.cat <- df.metric$category
df.norm <- apply(df.metric[, -"category"], 2, scale, 
                center = TRUE, scale = TRUE)
rownames(df.norm) <- df.metric$category

# filter for breakage types 
df.norm.filtered <- df.norm[which(grepl(
    pattern = "TFBS|DNA_Breakage|QM_parameters", x = gen.class)
),]

cor.out <- correlate.and.filter(dat = df.norm.filtered)
cor.out

fwrite(
    cor.out,
    paste0("../data/correlations/",
           "AllBreaks_QM_and_TFBS_", 
           kmer, "mer.csv")
)

################################################################################################################
# #' exploring correlation between QM parameters, biological and TFBS
# df.norm.filtered <- df.norm[which(grepl(
#     pattern = "TFBS|_ds_ss|WRN_loss", x = rownames(df.norm))
# ),]

# cor.out <- correlate.and.filter(dat = df.norm.filtered)
# cor.out

# fwrite(
#     cor.out,
#     paste0("../data/correlations/",
#            "QM_and_TFBS_kmer-", 
#            kmer, ".csv")
# )

cor.out <- cor.out %>% 
    tidyr::pivot_longer(
        cols = !tfbs, 
        names_to = "DNA_breakage", 
        values_to = "corr",
        values_drop_na = TRUE) %>% 
    dplyr::mutate(
        DNA_breakage = forcats::fct_inorder(DNA_breakage)
    )
     
# x-axis colourings
label.hex <- org_file$Category_Hex[match(
    levels(cor.out$DNA_breakage), org_file$Category_Main
)]
label.hex <- ifelse(is.na(label.hex), "#8F9CA7", label.hex)

p1 <- cor.out %>%
    dplyr::mutate(tfbs = stringr::str_remove(
        string = tfbs, pattern = "TFBS_")) %>% 
    dplyr::arrange(desc(tfbs)) %>% 
    dplyr::mutate(tfbs = forcats::fct_inorder(tfbs)) %>% 
    ggplot(aes(x = DNA_breakage, y = tfbs, fill = corr)) + 
    geom_tile(col = "grey") + 
    scale_fill_gradient2(
        low = "#2166ac",
        high = "#b2182b",
        midpoint = 0,
        limit = c(-1, 1),
        name = "Corr"
    ) + 
    theme_bw() + 
    theme_classic() + 
    suppressWarnings(theme(
        panel.grid.major.x = element_blank(),
        text = element_text(size = 20),
        # axis.text.x = element_blank()
        axis.text.x = element_text(
            angle = 45, 
            hjust = 1, 
            vjust = 1,
            colour = label.hex
        )
    )) + 
    labs(x = "", y = "")

ggsave(
    filename = paste0(
        "../figures/correlation_maps/", 
        "QM_and_TFBS_kmer-", 
        kmer, ".png"),
    plot = p1,
    height = 10,
    width = 9,
    dpi = 600
)

###########################################################################
# Extract sequence motifs from above filtered results
###########################################################################
unique.tfbs <- sort(unique(cor.out$tfbs))

# get the enriched k-mers for each tfbs
df.zscore <- df.metric[match(unique.tfbs, df.metric$category),]
df.zscore <- as.data.frame(df.zscore)
rownames(df.zscore) <- df.zscore$category
df.zscore <- subset(df.zscore, select = -c(category))

tfbs.dehof <- lapply(1:length(unique.tfbs), function(x){
    df.subset.t <- t(df.zscore[x,])
    arr.order <- order(df.subset.t, decreasing = TRUE)
    get.enriched <- which(df.subset.t[arr.order] >= 0)
    top.1.perc <- ceiling(length(df.subset.t)*0.01)
    if(length(get.enriched) > top.1.perc) get.enriched <- 1:top.1.perc

    # extract kmers
    enriched.kmers <- rownames(df.subset.t)[arr.order][get.enriched]

    # get delta Ehof 
    dehof.match <- df.dehof[match(enriched.kmers, df.dehof$kmer),]

    # add tfbs name
    dehof.match[, tfbs := gsub("TFBS_", "", unique.tfbs[x])]

    return(dehof.match)
})
tfbs.dehof <- rbindlist(tfbs.dehof)

tfbs.dehof

df.subset.t <- sort(t(df.zscore[1,]))
df.subset.t[which(df.subset.t >= 0)]

folder.names <- org_file$exp[match(unique.tfbs, org_file$Category_Main)]
files <- list.files(
    path = paste0("../data/kmertone/", folder.names),
    pattern = paste0("score_", kmer, "-mers.csv"),
    full.names = TRUE
)

lapply(1:length(files), function(x){
    dt <- fread(files[x], showProgress = FALSE)
    dt
})

# plot sequence logos

################################################################################################################
#' exploring correlation between cell free DNA and epigenome marks
df.norm.filtered <- df.norm[which(grepl(
    pattern = "TFBS|Epigenome|G4seq_map|Cell free DNA", x = labels)
),]

cor.out <- correlate.and.filter(
    dat = df.norm.filtered, 
    rows.filter = "TFBS|Epigenome|G4seq_map|seq"
)
cor.out

fwrite(
    cor.out,
    paste0("../data/correlations/",
           "CF-DNA_and_TFBS-Epigenome_kmer-", 
           kmer, ".csv")
)

cor.out <- cor.out %>% 
    tidyr::pivot_longer(
        cols = !tfbs, 
        names_to = "DNA_breakage", 
        values_to = "corr",
        values_drop_na = TRUE
    ) %>% 
    dplyr::mutate(
        DNA_breakage = forcats::fct_inorder(DNA_breakage)
    )
     
# x-axis colourings
label.hex <- org_file$Category_Hex[match(
    levels(cor.out$DNA_breakage), org_file$Category_Main
)]
label.hex <- ifelse(is.na(label.hex), "#8F9CA7", label.hex)

# cor.out <- 
# cor.out %>% 
#     dplyr::mutate(tfbs = forcats::fct_inorder(tfbs)) %>% 
#     dplyr::mutate(tfbs = stringr::str_remove_all(
#         string = tfbs,
#         pattern = "FAIREseq_|Epigenome_|Dnaseseq_|Chipseq_|ATACseq_"
#     )) %>% 
#     dplyr::mutate(DNA_breakage = stringr::str_remove_all(
#         string = DNA_breakage,
#         pattern = "Cancer_|Autoimmune_"),
#         DNA_breakage = stringr::str_to_title(DNA_breakage)
#     )

p1 <- cor.out %>% 
    ggplot(aes(x = DNA_breakage, y = tfbs, fill = corr)) + 
    geom_tile(col = "grey") + 
    scale_fill_gradient2(
        low = "#2166ac",
        high = "#b2182b",
        midpoint = 0,
        limit = c(-1, 1),
        name = "Corr"
    ) + 
    theme_bw() + 
    theme_classic() + 
    suppressWarnings(theme(
        text = element_text(size = 20),
        axis.text.x = element_text(
            angle = 45, 
            hjust = 1, 
            vjust = 1,
            colour = label.hex
        )
    )) + 
    labs(x = "", y = "")

ggsave(
    filename = paste0(
        "../figures/correlation_maps/", 
        "CF-DNA_and_TFBS-Epigenome_kmer-", 
        kmer, ".png"),
    plot = p1,
    height = 8,
    width = 11,
    dpi = 600
)