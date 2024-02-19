suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(gridExtra)))
suppressPackageStartupMessages(suppressWarnings(library(cluster)))
suppressPackageStartupMessages(suppressWarnings(library(mclust)))
suppressPackageStartupMessages(suppressWarnings(library(gplots)))
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

kmer=8

# get rmsd tract for kmer
df <- fread(paste0("../data/ranges/kmer_", kmer, 
            "_Ranges_cutoffs_from_clustering_all-exp.csv"))

# add breakage category
org.file <- fread("../../data/org_file.csv")
org.file[, exp := paste0(`Fragmentation type`, "/", `Experiment folder`)]
inds <- match(df$exp, org.file$exp)
category.labels <- org.file$Category_Main[inds]
category.general.labels <- org.file$Category_general[inds]
category.general.hex <- org.file$Category_Hex[inds]
df[, category := category.labels]
df[, `Breakage Type` := category.general.labels]
df[, hex.cols := category.general.hex]

########################################################################
# Average sequence context influence per curve
########################################################################
# for manual plot colourings
if(is.null(df$`long.range-2`)) df <- dplyr::mutate(df, `long.range-2` = NA)
if(is.null(df$`mid.range-2`)) df <- dplyr::mutate(df, `mid.range-2` = NA)
if(is.null(df$`short.range-2`)) df <- dplyr::mutate(df, `short.range-2` = NA)

df.mean.ranges <- as_tibble(df) %>% 
    dplyr::group_by(exp, rowid) %>% 
    dplyr::mutate(across(
            dplyr::matches("long|mid|short"),
            ~ tidyr::replace_na(., 0)
        )
    ) %>% 
    dplyr::mutate(
        long.range = dplyr::case_when(
            rowid == "ranges" ~ max(c(`long.range-1`, `long.range-2`), na.rm = TRUE),
            TRUE ~ sum(c(`long.range-1`, `long.range-2`), na.rm = TRUE)
        ),
        mid.range = dplyr::case_when(
            rowid == "ranges" ~ max(c(`mid.range-1`, `mid.range-2`), na.rm = TRUE),
            TRUE ~ sum(c(`mid.range-1`, `mid.range-2`), na.rm = TRUE)
        ),
        short.range = dplyr::case_when(
            rowid == "ranges" ~ max(c(`short.range-1`, `short.range-2`), na.rm = TRUE),
            TRUE ~ sum(c(`short.range-1`, `short.range-2`), na.rm = TRUE)
        )
    ) %>% 
    mutate(
        short.range = dplyr::case_when(
            rowid == "contribution" ~ short.range*100,
            TRUE ~ short.range
        )
    ) %>%
    dplyr::select(-c(
        dplyr::contains(".range-1"),
        dplyr::contains(".range-2"))
    ) %>%  
    dplyr::group_by(`Breakage Type`, category, hex.cols, rowid) %>% 
    dplyr::summarise(
        mean.lr = mean(long.range, na.rm = TRUE),
        mean.mr = mean(mid.range, na.rm = TRUE),
        mean.sr = mean(short.range, na.rm = TRUE),
        .groups = "keep"
    ) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(`Breakage Type` = as.factor(`Breakage Type`)) 

filter.df.onlyranges <- function(cols, col.label){
    label <- switch(col.label,
        "sr" = "Short range",
        "mr" = "Medium range",
        "lr" = "Long range"
    )

    df.for.plotting <- df.mean.ranges %>% 
        dplyr::filter(rowid == "ranges") %>%
        dplyr::select(category, `Breakage Type`, hex.cols, {{cols}}) %>% 
        tidyr::drop_na() %>% 
        dplyr::arrange({{cols}}, desc(category)) %>% 
        dplyr::mutate(category = forcats::fct_inorder(category)) %>% 
        dplyr::filter({{cols}} > 0)

    p1 <- df.for.plotting %>% 
        ggplot(aes(x = {{cols}}, y = category)) + 
        geom_segment(aes(xend = 0, yend = category), 
            col = "darkgrey",
            linewidth = 0.6) +
        geom_point(aes(x = {{cols}}, y = category),
            color = df.for.plotting$hex.cols, 
            size = 4) + 
        theme_bw() +
        theme_classic() + 
        theme(axis.text.x = element_text(size = 15)) +
        labs(
            x = paste0(label, " influence (avg)"),
            y = ""
        )
        
    return(p1)
}

if(kmer == 8){
    p1 <- filter.df.onlyranges(cols = mean.sr, col.label = "sr")
    p2 <- filter.df.onlyranges(cols = mean.mr, col.label = "mr")
    p3 <- filter.df.onlyranges(cols = mean.lr, col.label = "lr")
    pdf(
        paste0("../figures/seq_influence/", kmer, "-mer-", 
            "avg-ranges-ranked.pdf"),
        height = 18, width = 9
    )
    gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
    plot.saved <- dev.off()
}

########################################################################
# Average sequence context influence per curve vs. AUC or gaussian peaks
########################################################################
df.for.plotting <- df.mean.ranges %>% 
    dplyr::filter(rowid == "peak.intensity" | rowid == "ranges") %>% 
    dplyr::mutate(rowid = ifelse(rowid == "peak.intensity", "contribution", "ranges")) %>% 
    dplyr::select(category, `Breakage Type`, rowid, hex.cols, dplyr::starts_with("mean")) %>% 
    dplyr::mutate(category = forcats::fct_inorder(category)) %>% 
    tidyr::pivot_longer(
        cols = -c(category, `Breakage Type`, hex.cols, rowid),
        names_to = "key",
        values_to = "value"
    ) %>% 
    dplyr::filter(value != 0) %>% 
    tidyr::pivot_wider(
        names_from = rowid,
        values_from = value
    )

df.for.plotting$key <- ifelse(
    grepl(pattern = "sr", x = df.for.plotting$key), "Short range",
    ifelse(grepl(pattern = "mr", x = df.for.plotting$key), "Medium range",
    ifelse(grepl(pattern = "lr", x = df.for.plotting$key), "Long range", 
    df.for.plotting$key
)))

# contribution plot
p1 <- df.for.plotting %>%
    dplyr::mutate(key = factor(key, 
        levels = c("Short range", "Medium range", "Long range"))) %>% 
    ggplot(aes(x = ranges, y = contribution)) + 
    geom_point(col = df.for.plotting$hex.cols, size = 5, alpha = 1) +        
    facet_wrap(vars(key), nrow = 1, scales = "free_x") + 
    theme_bw() + 
    labs(x = "", y = "") +   
    coord_cartesian(
        xlim = c(0, NA), 
        ylim = c(0, 100)
    ) +
    theme(
        text = element_text(size = 25),
        plot.margin = grid::unit(c(1,1,1,-1), "cm")
    ) 
ggsave(
    filename = paste0("../figures/seq_influence/", kmer, "-mer-", 
            "avg-ranges-ranked_with-GaussianPeaksStats.pdf"),
    plot = p1,
    height = 7, width = 13
)

# contribution plot split into short, mid and long range
df.for.plotting <- df.for.plotting %>%
    dplyr::mutate(key = factor(key, 
        levels = c("Short range", "Medium range", "Long range"))
    )

df.for.plotting.ranges <- df.for.plotting %>% 
    dplyr::group_by(key) %>% 
    dplyr::summarise(
        min.range = min(ranges, na.rm = TRUE),
        max.range = max(ranges, na.rm = TRUE)
    )
df.for.plotting.SR <- df.for.plotting %>% dplyr::filter(key == "Short range")
df.for.plotting.MR <- df.for.plotting %>% dplyr::filter(key == "Medium range")
df.for.plotting.LR <- df.for.plotting %>% dplyr::filter(key == "Long range")
df.for.plotting.ranges[df.for.plotting.ranges$key == "Short range","min.range"] <- 0

#' @description 
#' Plot each range of sequence influence manually
#' @return ggplot object.
plot.each.range <- function(dat, x.range = c(0, 35), y.ticks = FALSE){
    # manual labelling
    p1 <- dat %>% 
        ggplot(aes(x = ranges, y = contribution)) + 
        geom_point(col = dat$hex.cols, size = 5, alpha = 0.8) +        
        facet_wrap(vars(key), nrow = 1, scales = "free_x") + 
        theme_bw() + 
        labs(x = "", y = "") +   
        coord_cartesian(
            xlim = x.range, 
            ylim = c(0, 100)
        ) +
        theme(
            text = element_text(size = 30),
            plot.margin = grid::unit(c(1,1,1,0), "cm"),
            strip.text = element_blank(),
            strip.background = element_rect(fill = "white", colour = "white")
        )
    
    if(!y.ticks) p1 <- p1 + theme(axis.text.y = element_blank())

    return(p1)
}

# use ranges here as guidance
df.for.plotting.ranges

p1 <- plot.each.range(dat = df.for.plotting.SR, x.range = c(0, 22))
p2 <- plot.each.range(dat = df.for.plotting.MR, x.range = c(25, 130))
p3 <- plot.each.range(dat = df.for.plotting.LR, x.range = c(140, 2000))

ggsave(
    filename = paste0("../figures/seq_influence/", kmer, "-mer-", "demo_3.png"),
    plot = p1,
    height = 8, width = 8
)

png(
    file = paste0("../figures/seq_influence/", kmer, "-mer-", 
                  "avg-ranges-ranked_with-GaussianPeaksStats_v2.png"),
    height = 7, width = 20,
    units = "in", res = 600 
)
gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
plot.saved <- dev.off()

pdf(
    file = paste0("../figures/seq_influence/", kmer, "-mer-", 
                  "avg-ranges-ranked_with-GaussianPeaksStats_v2.pdf"),
    height = 7, width = 20
)
gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
plot.saved <- dev.off()

# with ticks
p1 <- plot.each.range(dat = df.for.plotting.SR, x.range = c(0, 22), y.ticks = TRUE)
p2 <- plot.each.range(dat = df.for.plotting.MR, x.range = c(25, 130))
p3 <- plot.each.range(dat = df.for.plotting.LR, x.range = c(140, 2000))

png(
    file = paste0("../figures/seq_influence/", kmer, "-mer-", 
                  "avg-ranges-ranked_with-GaussianPeaksStats_v2_ytext.png"),
    height = 7, width = 20,
    units = "in", res = 600 
)
gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
plot.saved <- dev.off()

# get maximum value for each range effect by type
get_max_val_by_type <- function(dat, which_range){
    max_short_by_type <- dat %>% 
        dplyr::group_by(`Breakage Type`) %>% 
        dplyr::summarise(max_range = max(ranges, na.rm = TRUE)) %>% 
        dplyr::rename_with(~c("break_type", "max_range")) %>% 
        dplyr::mutate(
            break_type = gsub(" ", "_", break_type),
            break_type = dplyr::case_when(
                grepl("Cell_free|Mechanical|Natural", break_type) ~ "High_frequency",
                TRUE ~ break_type
            )
        ) %>% 
        dplyr::group_by(break_type) %>% 
        dplyr::summarise(max_range = max(max_range, na.rm = TRUE)) %>% 
        dplyr::arrange(break_type) %>% 
        dplyr::rename_with(~c("break_type", which_range))
    return(max_short_by_type)
}

df_max_range_by_type <- dplyr::left_join(
    x = get_max_val_by_type(
        dat = df.for.plotting.SR, which_range = "short"
    ),
    y = get_max_val_by_type(
        dat = df.for.plotting.MR, which_range = "medium"
    ),
    by = "break_type",
    keep = FALSE
    ) %>% 
    dplyr::left_join(.,     
        y = get_max_val_by_type(
            dat = df.for.plotting.LR, which_range = "long"
        ),
        by = "break_type",
        keep = FALSE
    ) %>% 
    as.data.table()

max_range_effects <- df.for.plotting %>% 
    dplyr::group_by(`Breakage Type`, key) %>% 
    dplyr::summarise(max_range = max(ranges, na.rm = TRUE)) %>% 
    dplyr::rename_with(~c("Breakage_type", "Range_effect", "Max_range")) %>% 
    dplyr::mutate(Range_effect = stringr::str_remove(
        string = Range_effect, pattern = " range")
    ) %>% 
    suppressMessages()

if(kmer == 8){
    fwrite(
        df_max_range_by_type,
        file = "../data/ranges/MaxValuesFromClustersByType.csv"
    )
    
    fwrite(
        df_max_range_by_type,
        file = paste0(
            "../../../04_DNAFragility/data/",
            "range_effects/MaxValuesFromClustersByType.csv"
        )
    )

    fwrite(
        max_range_effects,
        file = paste0(
            "../../../04_DNAFragility/data/",
            "range_effects/MaxRangeEffects.csv"
        )
    )
}

########################################################################
# Kmeans clustering between points on linear combination Gaussian curve
########################################################################
# import gaussian curves
files <- list.files(
    path = "../data",
    pattern = paste0("gaussian_curve_vals_linear_combinations_kmer_", kmer, ".csv"),
    full.names = TRUE,
    recursive = TRUE
)
files <- files[!grepl(pattern = "AverageRMSDFromClusters", x = files)]
file.names <- stringr::str_remove_all(
    string = files,
    pattern = paste0("../data/|/gaussian_curve_vals_linear_combinations_kmer_", kmer, ".csv")
)

# renamings
col.names <- stringr::str_remove_all(
    string = unique(file.names),
    pattern = "00_Breakage/"
)
inds <- match(unique(file.names), org.file$exp)
category.labels <- org.file$Category_Main[inds]
category.general.labels <- org.file$Category_general[inds]
category.general.hex <- org.file$Category_Hex[inds]

df.to.import <- data.table(
    exp = unique(file.names),
    exp.clean = col.names,
    category = category.labels,
    breakage.type = category.general.labels,
    hex.cols = category.general.hex,
    file = files
)
to.import <- unname(split(df.to.import, df.to.import$category))

res <- lapply(1:length(to.import), function(x){
    dfs <- lapply(1:nrow(to.import[[x]]), function(i){
        temp <- fread(to.import[[x]]$file[i], showProgress = FALSE)
        setnames(temp, to.import[[x]]$exp[i])
        return(temp)
    })        
    dfs <- do.call(cbind, dfs)
    dfs <- as.data.table(rowMeans(dfs, na.rm = TRUE))
    setnames(dfs, to.import[[x]]$category[1])
    return(dfs)
})
res <- do.call(cbind, res)
res.original <- t(res)

perform_hclust <- function(res_mat, rm_labels = FALSE){
    if(rm_labels){        
        res_mat <- res_mat[!grepl("enzyme", rownames(res_mat)),]
    }

    # normalise and rename
    res <- apply(res_mat, 2, scale, 
                center = TRUE, scale = TRUE)
    rownames(res) <- rownames(res_mat)

    # euclidean distances
    res.dist <- dist(res) %>% suppressWarnings()

    # finding optimal number of clusters for hierarchical clustering
    suppressPackageStartupMessages(library(factoextra))

    max.ks <- 15
    p1 <- fviz_nbclust(res, hcut, method = "wss", 
                    k.max = max.ks, linecolor = "#6D7997") + 
        theme(text = element_text(size = 15)) + 
        coord_cartesian(xlim = c(1,max.ks))
    p2 <- fviz_nbclust(res, hcut, method = "silhouette", 
                    k.max = max.ks, linecolor = "#6D7997") + 
        coord_cartesian(xlim = c(1,max.ks)) + 
        theme(text = element_text(size = 15)) + 
        labs(title = "", x = "")

    png(paste0(
        "../figures/seq_influence/", kmer, 
        "-mer-allGaussianCurves_kmeans_eval.png"),
        height = 8, width = 12,
        units = "in", res = 600 
    )
    gridExtra::grid.arrange(p1, p2, ncol = 2)
    plot.saved <- dev.off()
    
    # hierarchical clustering
    suppressPackageStartupMessages(library(dendextend))

    max.ks <- 6
    res.hclust <- hclust(res.dist, method = "ward.D2")
    res.dendo <- as.dendrogram(res.hclust)
    res.clusters <- cutree(res.hclust, k = max.ks)

    # colour by breakage type
    inds <- match(rownames(res), org.file$Category_Main)
    category.labels <- org.file$Category_Main[inds]
    category.general.labels <- org.file$Category_general[inds]
    category.general.hex <- org.file$Category_Hex[inds]
    col.df <- data.table(
        exp = rownames(res),
        category = category.labels,
        breakage.type = category.general.labels,
        hex.cols = category.general.hex
    )
    col.df[, clusters := unname(res.clusters)[match(
        col.df$exp, names(res.clusters))]]

    # order by dendogram
    col.df <- col.df[order(match(
        exp, col.df$exp[order.dendrogram(res.dendo)]
    ))]
    labels_colors(res.dendo) <- col.df$hex.cols

    # # colour cluster
    # fac <- factor(rep_len(c("#cecece", "#9d9d9d"), length.out = max.ks))
    # num <- as.numeric(fac)
    # new.num <- num[unname(res.clusters)]

    # plot
    pdf(
        paste0("../figures/seq_influence/", kmer, 
            "-mer-allGaussianCurves_hclust_zoomed_outer.pdf"), 
        height = 10, width = 7
    )
    par(mar = c(4,1,1,22))
    plot(res.dendo %>% set("labels_cex", c(1.2, 1.2)), horiz = TRUE, xlim = c(350, 0))
    plot.save <- dev.off()

    pdf(
        paste0("../figures/seq_influence/", kmer, 
            "-mer-allGaussianCurves_hclust_zoomed.pdf"), 
        height = 10, width = 7
        # units = "in", res = 600 
    )
    par(mar = c(4,1,1,22))
    plot(res.dendo %>% set("labels_cex", c(1.2, 1.2)), horiz = TRUE, xlim = c(130, 0))
    plot.save <- dev.off()

    pdf(
        paste0("../figures/seq_influence/", kmer, 
            "-mer-allGaussianCurves_hclust_withRECT_zoomed.pdf"), 
        height = 10, 
        width = 7
    )
    par(mar = c(4,1,1,22))
    plot(res.dendo, horiz = TRUE)
    rect.dendrogram(res.dendo, k=max.ks,horiz = TRUE, xlim = c(50, 0))
    plot.save <- dev.off()

    res.clusters <- setNames(col.df$clusters, col.df$exp)
    return(res.clusters)
}

#' removing enzymes and the 'elbow method' tells us 3 clusters,
#' and then the enzymes all get their own cluster. 
res.clusters <- perform_hclust(res_mat = res.original, rm_labels = FALSE)
# perform_hclust(res_mat = res.original, rm_labels = TRUE)

########################################################################
# RMSD plots + Gaussian curve fittings based on average curves 
# obtained above from kmeans clustering
########################################################################
if(kmer == 8){
    # assign names with clusters
    res.cluster.df <- tibble(
            exp =  names(res.clusters),
            cluster = unname(res.clusters)
        ) %>% 
        dplyr::mutate(cluster = match(cluster, unique(cluster))) %>% 
        as.data.table()
    setorder(res.cluster.df, cluster)
    unique.clusters <- unique(res.cluster.df$cluster)

    #' note to self: make sure you check what, for example, which exp were in 
    #' cluster 1 as per the processing below and where you put that plot on the 
    #' above hierarchical clustering plot.
    # res.cluster.df

    source("../lib/FitCurves.R")
    fitted_curves <- lapply(rev(unique(res.cluster.df$cluster)), function(x){
        # for each cluster, get the average RMSD vector
        exps <- res.cluster.df[cluster == x]$exp
        row.ind <- match(exps, rownames(res.original))

        if(length(row.ind) > 1){
            avg.rmsd <- colMeans(res.original[row.ind, ], na.rm = TRUE)
        } else {
            avg.rmsd <- res.original[row.ind, ]
        }

        # get RMSD plot and fit Gaussian curves
        set.seed(1234)
        fit_curves <- FitCurves$new(
            chr = 1,
            rmsd_values = list("chr1" = avg.rmsd),
            return_vals = TRUE
        )
        fit_curves$generate_rmsd_plots(k = kmer, annot_plots = FALSE, per_chr = TRUE)

        avg.rmsd.df <- rbind(
            data.table(vals = avg.rmsd, col.label = "RMSD"),
            fit_curves$results$chr1$curve_vals,
            data.table(
                vals = fit_curves$results$chr1$all_gc$V1, 
                col.label = "curve.all"
            )
        )
        avg.rmsd.df[, cluster := x]

        return(list(
            avg.rmsd.df,
            fit_curves$results$chr1$df %>% 
                dplyr::mutate(cluster = x)
        ))
    })

    avg.rmsd.df <- rbindlist(sapply(fitted_curves, `[`, 1))
    avg.rmsd.stats <- do.call(rbind, sapply(fitted_curves, `[`, 2))

    # save result
    dir.create(path = "../data/AverageRMSDFromClusters", showWarnings = FALSE)

    saveRDS(
        avg.rmsd.df, 
        file = "../data/AverageRMSDFromClusters/Cluster_RMSD.Rdata"
    )
    fwrite(
        avg.rmsd.stats, 
        file = "../data/AverageRMSDFromClusters/Cluster_RMSD_stats.csv"
    )

    hash.map <- df.for.plotting.ranges %>% 
        dplyr::mutate(Ranges = case_when(
            max.range == min(max.range, na.rm = TRUE) ~ "short.range",
            max.range == median(max.range, na.rm = TRUE) ~ "mid.range",
            max.range == max(max.range, na.rm = TRUE) ~ "long.range")) %>% 
        dplyr::pull(max.range, Ranges)

    res.clusters <- lapply(unique(avg.rmsd.df$cluster), function(x){
        df <- avg.rmsd.df[cluster == x]
        df_stats <- as.data.table(avg.rmsd.stats[avg.rmsd.stats$cluster == x,])

        df <- as_tibble(df) %>% 
            dplyr::select(-cluster) %>% 
            split(f = as.factor(.$col.label)) %>% 
            lapply(., function(x){
                x %>% 
                    dplyr::select(-col.label) %>%
                    dplyr::rename_with(~unique(x$col.label))
            }) %>% 
            do.call(dplyr::bind_cols, .)

        # if only two ranges present, then add longer range influence as NA
        all_cols <- which(grepl(pattern = "one|two|three", x = colnames(df)))
        if(length(all_cols) < 3){
            df <- df %>% dplyr::mutate(curve.three = NA_integer_)
        } 

        # assign curves 1-3 into short, medium or long range
        assigning.ranges <- df_stats[rowid == "ranges"]
        cols.of.int <- which(grepl(pattern = "curve", x = colnames(assigning.ranges)))
        assigning.ranges <- assigning.ranges[, ..cols.of.int]

        range_assignments <- as_tibble(assigning.ranges) %>% 
            tidyr::pivot_longer(
                cols = everything(),
                names_to = "key",
                values_to = "value"
            ) %>% 
            dplyr::mutate(
                Cluster = case_when(
                    value <= min(hash.map) ~ "Short range",
                    value >= min(hash.map) & value <= median(hash.map) ~ "Medium range",
                    value > median(hash.map) & value <= max(hash.map) ~ "Longer range",
                    TRUE ~ NA_character_
                )
            )

        # if any Cluster has multiple occurrences, then take max
        if(any(duplicated(range_assignments$Cluster))){
            to_merge <- range_assignments %>% 
                dplyr::group_by(Cluster) %>% 
                dplyr::mutate(value = max(value))

            which_dup <- which(
                duplicated(to_merge$Cluster) | duplicated(to_merge$Cluster, fromLast = TRUE)
            )
            which_dup <- match(c("rowid", to_merge$key[which_dup]), colnames(df_stats))
            
            df_stats <- as_tibble(df_stats[, ..which_dup]) %>% 
                dplyr::mutate(cluster = x) %>% 
                dplyr::mutate(
                    to_rename = case_when(
                        rowid %in% c("ranges", "SD") ~ dplyr::select(., starts_with("curve")) %>% 
                            purrr::pmap_dbl(~ max(c(...), na.rm = TRUE)),
                        TRUE ~ select(., starts_with("curve")) %>% 
                            rowSums(na.rm = TRUE)
                    )
                ) %>% 
                dplyr::select(rowid, to_rename, cluster) %>% 
                dplyr::rename_with(~c(
                    "rowid", 
                    to_merge$key[duplicated(to_merge$Cluster, fromLast = TRUE)],
                    "cluster"
                ))

            col_exists <- which(grepl(pattern = "curve.two", x = colnames(df_stats)))
            if(length(col_exists) == 0){
                df_stats <- dplyr::mutate(df_stats, curve.two = NA_real_)

                to_merge$key[match("curve.two", to_merge$key)] <- "curve.two"
                to_merge$value[match("curve.two", to_merge$key)] <- NA_real_
                to_merge$Cluster[match("curve.two", to_merge$key)] <- NA_real_
            }

            col_exists <- which(grepl(pattern = "curve.three", x = colnames(df_stats)))
            if(length(col_exists) == 0){
                df_stats <- dplyr::mutate(df_stats, curve.three = NA_real_)

                to_merge$key[match("curve.three", to_merge$key)] <- "curve.three"
                to_merge$value[match("curve.three", to_merge$key)] <- NA_real_
                to_merge$Cluster[match("curve.three", to_merge$key)] <- NA_real_
            }

            df_stats <- df_stats %>% 
                dplyr::select(rowid, dplyr::starts_with('curve'), cluster)

            range_assignments <- to_merge 
        }
        df_stats <- as.data.table(df_stats)

        unused_range <- setdiff(
            c("Short range", "Medium range", "Longer range"), 
            unique(range_assignments$Cluster)
        )

        range_assignments <- range_assignments %>% 
            dplyr::mutate(
                Cluster = ifelse(is.na(Cluster), unused_range, Cluster)
            ) %>% 
            dplyr::pull(Cluster, key)

        # some DSB processes may have 2 separate sequence-context effects within the same
        # categorised range window. In that case, re-label the longest one as "Longer range"
        # and point out in manuscript the difference
        all_cols <- which(grepl(pattern = "Medium", x = range_assignments))
        if(length(all_cols) > 1){
            range_assignments[3] <- "Longer range"
        }

        # rename curves
        all_cols <- which(grepl(pattern = "one|two|three", x = colnames(df)))
        col_rename <- as.character(range_assignments[colnames(df)[all_cols]])
        colnames(df)[all_cols] <- col_rename
        df <- df %>% dplyr::mutate(cluster = x) %>% as.data.table()
        setcolorder(df, c("curve.all", "Short range", "Medium range", "Longer range", "cluster"))

        limits <- nrow(df)/2-1
        df$x <- (-limits):(nrow(df)-limits-1)

        # repeat for stats table
        all_cols <- which(grepl(pattern = "one|two|three", x = colnames(df_stats)))
        col_rename <- as.character(range_assignments[colnames(df_stats)[all_cols]])
        colnames(df_stats)[all_cols] <- col_rename

        df_stats <- df_stats[rowid == "peak.intensity" | rowid == "ranges"]
        setorder(df_stats, rowid)
        df_stats[, rowid := factor(
            c("Contribution to full range, %", "Range of influence"), 
            levels = c("Range of influence", "Contribution to full range, %"))]
        setcolorder(df_stats, c("rowid", "Short range", "Medium range", "Longer range", "cluster"))
        return(list(df, df_stats))
    })
    df_rmsd <- rbindlist(sapply(res.clusters, `[`, 1))
    df_stats <- rbindlist(sapply(res.clusters, `[`, 2))

    # for plotting the RMSD curves
    df_rmsd <- as_tibble(df_rmsd) %>% 
        tidyr::pivot_longer(
            cols = -c(x, cluster),
            names_to = "Key",
            values_to = "Value"
        ) %>% 
        dplyr::mutate(
            Key = ifelse(Key == "curve.all", "Full range", Key),
            hex = dplyr::case_when(
                Key == "RMSD" ~ "darkgrey",
                Key == "Short range" ~ "#619B61",
                Key == "Medium range" ~ "#355c7d",
                Key == "Longer range" ~ "#BB357E",
                Key == "Full range" ~ "#f8b195"
            ),
            alpha = dplyr::case_when(
                Key == "RMSD" ~ 0.5,
                TRUE ~ 1
            ),
            lw = dplyr::case_when(
                Key == "RMSD" ~ 0.5,
                Key == "Short range" ~ 1,
                Key == "Medium range" ~ 1,
                Key == "Longer range" ~ 1,
                Key == "Full range" ~ 1.2
            ),
            Key = factor(Key, levels = c(
                "RMSD",
                "Short range", 
                "Medium range", 
                "Longer range", 
                "Full range"
                )
            ),
            cluster = paste0("Cluster: ", cluster)
        ) %>% 
        tidyr::drop_na()

    # Creating a mapping of old levels to new titles
    cluster_titles <- setNames(
        paste("Cluster", seq_along(rev(unique(df_rmsd$cluster)))), 
        unique(df_rmsd$cluster)
    )

    p1 <- df_rmsd %>% 
        dplyr::mutate(
            cluster = factor(cluster, levels = names(cluster_titles)),
            cluster = factor(cluster_titles[as.character(cluster)])
        ) %>% 
        ggplot(aes(x = x, y = Value, col = hex, alpha = alpha)) + 
        geom_line(linewidth = 1.5) +
        facet_wrap(vars(cluster), ncol = 1, scales = "free_y") + 
        scale_color_identity() + 
        theme_bw() + 
        theme_classic() + 
        coord_cartesian(xlim = c(-100, 100)) + 
        theme(
            text = element_text(size = 20),
            legend.position = "none"
            # strip.text = element_blank()
        ) + 
        labs(
            x = "Position away from breakpoint, bp",
            y = "RMSD"
        )

    dir.create(path = "../figures/AverageRMSDFromClusters", showWarnings = FALSE)
    ggsave(
        filename = paste0(
            "../figures/AverageRMSDFromClusters/",
            "All_curve_fits-kmer_", kmer, ".pdf"
        ),
        plot = p1,
        height = 11, width = 7
    )

    # line plot for short, Medium and long range statistics
    res.clusters <- as_tibble(df_stats) %>% 
        tidyr::pivot_longer(
            !c(rowid, cluster), 
            names_to = "range", 
            values_to = "values"
        ) %>% 
        dplyr::mutate(hex.col = 
            dplyr::case_when(
                range == "Short range" ~ "#619B61",
                range == "Medium range" ~ "#355c7d",
                range == "Longer range" ~ "#BB357E",
                TRUE ~ NA_character_
            )
        ) %>% 
        dplyr::mutate(values = ifelse(is.na(values), 0, values))

    # some DSB processes may have 2 separate sequence-context effects within the same
    # categorised range window. In that case, re-label the longest one as "Longer range"
    # and point out in manuscript the difference. Here, we know that cluster 2 has 
    # two medium range and no categorised long-range effect.
    # res.clusters[(res.clusters$cluster == 4) & (res.clusters$range == "Longer range"), "hex.col"] <- "#000000"

    res.clusters.left <- res.clusters %>%
        dplyr::filter(rowid == "Range of influence") %>% 
        dplyr::mutate(values = ifelse(values == 0, NA_integer_, values))
        
    res.clusters.right <- res.clusters %>%
        dplyr::filter(rowid == "Contribution to full range, %")

    # range of influence - left plot
    p1.left <- res.clusters.left %>% 
        ggplot(aes(x = values, y = cluster)) + 
        geom_segment(aes(
            xend = 0, yend = cluster), 
            col = "darkgrey",
            alpha = 0.8) +
        geom_point(col = res.clusters.left$hex.col, size = 5, alpha = 1) + 
        coord_cartesian(
            ylim = c(0.8, length(unique(res.cluster.df$cluster))+0.2)
        ) +
        facet_wrap(vars(rowid), ncol = 4, scales = "free_x") + 
        theme_bw() + 
        theme(
            text = element_text(size = 20),
            axis.text.y = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")
        ) +
        scale_y_continuous(breaks = 0:(length(unique(res.cluster.df$cluster))+1)) +
        labs(x = "", y = "")


    # contribution to full range - right plot
    p1.right <- res.clusters.right %>% 
        ggplot(aes(x = cluster, y = values, label = signif(values, 2))) + 
        geom_bar(position = "stack", stat = "identity", fill = res.clusters.right$hex.col) + 
        coord_flip(xlim = c(0.8, length(unique(res.cluster.df$cluster))+0.2)) + 
        facet_wrap(vars(rowid), ncol = 4, scales = "free_x") + 
        theme_bw() + 
        theme(
            text = element_text(size = 20),
            axis.text.y = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")
        ) +
        scale_x_continuous(
            breaks = 0:(length(unique(res.cluster.df$cluster))+1)
        ) +
        labs(x = "", y = "")

    temp = df_rmsd %>% 
        dplyr::mutate(
            cluster = factor(cluster, levels = c(
                rev(unique(df_rmsd$cluster))
            ))
        )

    p1 <- p1 + 
        theme(axis.text.y = element_blank()) +
        labs(x = "", y = "")

    pdf(
        file = paste0(
            "../figures/AverageRMSDFromClusters/",
            "All_curve_fits_summarystats-kmer_", kmer, ".pdf"
        ), 
        height = 10, width = 11
        # units = "in", res = 600 
    )
    gridExtra::grid.arrange(p1, p1.left, p1.right, ncol = 3)
    plot.saved <- dev.off()

    res.clusters.right %>% 
        dplyr::filter(range == "Short range")
}