# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
k <- as.integer(args[2])

# my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/01_RMSD/scripts/"
# k=4
setwd(my.path)

# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

# import RMSD data
all.files <- list.files(
  path = "../data/", 
  pattern = paste0("key_stats_kmer_", k, " *"), 
  recursive = TRUE 
)

results <- lapply(all.files, function(file){
  out <- fread(paste0("../data/", file))
  out <- as_tibble(out)

  if(!("curve.three" %in% colnames(out))){
    out <- cbind(out, rep(NA_real_, 4))
    colnames(out)[length(colnames(out))] <- "curve.three"
  }

  return(out)
}) 

results <- do.call(rbind, results)

# Flatten all ranges; plot histogram/density plot
df <- results %>% 
    dplyr::filter(rowid == "ranges") %>% 
    dplyr::select(exp, curve.one, curve.two, curve.three) %>% 
    tidyr::gather(-exp, key = "Curve", value = "Value") %>% 
    dplyr::filter(!is.na(Value))

# df %>% 
#     ggplot(aes(x = Value)) +
#     geom_histogram(
#         aes(y = ..density..),
#         fill = "skyblue",
#         colour = "black",
#         bins = 90,
#         alpha = 0.4
#     ) + 
#     scale_x_log10() + 
#     geom_density() + 
#     labs(
#         x = "Ranges",
#         y = "Density"
#     )

# Try to fit binomial curve to this
Clustering <- function(dat, log_scale){
    dat$log_Value <- log10(dat$Value)

    # Clustering
    if(log_scale){
        dat.dist <- dist(dat$log_Value)
    } else {
        dat.dist <- dist(dat$Value)
    }

    clusts <- hclust(dat.dist, method = "complete")
    cl_members <- cutree(clusts, k = 3)

    dat <- dat %>% 
        mutate(Cluster = as.numeric(cl_members))
        
    hash_map <- dat %>% 
        group_by(Cluster) %>% 
        summarise(
            lower.end = min(Value),
            upper.end = max(Value)
        ) %>% 
        mutate(Ranges = case_when(
            upper.end == min(upper.end) ~ "short.range",
            upper.end == median(upper.end) ~ "mid.range",
            upper.end == max(upper.end) ~ "long.range",
        )) %>% 
        pull(Ranges, Cluster)

    dat <- dat %>% 
        mutate(Cluster = case_when(
           Cluster == names(hash_map)[1] ~ unname(hash_map)[1],
           Cluster == names(hash_map)[2] ~ unname(hash_map)[2],
           Cluster == names(hash_map)[3] ~ unname(hash_map)[3]
        ))
    
    if(log_scale) df <<- dat 

    dat.plot <- dat %>% 
        ggplot(aes(x = Value)) +
        geom_histogram(aes(
            y = ..density..,
            fill = Cluster, 
            group = Cluster), 
            colour = "black",
            bins = 130
        ) +
        # scale_fill_discrete(labels = label) +   
        labs(
            title = paste("Sequence influences clustered into", length(unique(dat$Cluster)), "separate ranges"),
            subtitle = ifelse(log_scale, "Clustering based on log-scaled values", "Clustering based on true values"),
            x = "Ranges (raw values)",
            y = "Density in Cluster"
        )
        
    ggsave(
        plot = dat.plot,
        filename = paste0("../figures/kmer_", k, ifelse(log_scale, "_LOG_", "_"), "flattened_ranges_histogram.pdf"),
        width = 10, height = 7
    )

    if(!log_scale){
        # Hierarchical clustering plot
        pdf(paste0("../figures/kmer_", k, ifelse(log_scale, "_LOG_", "_"), "flattened_ranges_clustering.pdf"))
        plot(
            clusts, 
            hang = -1, 
            cex = 0.2,
            main = "Hierarchical clustering on flattened ranges of sequence influences",

            xlab = "Ranges", 
            sub = NA
        )

        rect.hclust(
            clusts, 
            k = max(cl_members), 
            which = seq_len(max(cl_members)), 
            border = seq_len(max(cl_members)), 
            cluster = cl_members
        )
        save.plot <- dev.off()
    }
}

Clustering(dat = df, log_scale = TRUE)
Clustering(dat = df, log_scale = FALSE)
 
# find cut-off values for each range
cutoff.ranges <- df %>% 
    group_by(Cluster) %>% 
    summarise(
        lower.end = min(Value),
        upper.end = max(Value)
    )

fwrite(
    x = cutoff.ranges, 
    file = paste0("../data/kmer_", k, "_Ranges_cutoffs_from_clustering.csv")
)

files.to.import <- list.files(
    path = "../data/", 
    pattern = "_Ranges_cutoffs"
)

kmer.val <- as.integer(str_extract(
    string = files.to.import, 
    pattern = "[:digit:]"
))

if(k == 8){
    out <- lapply(1:length(files.to.import), function(i){
        out <- fread(paste0("../data/", files.to.import[i]))
        out <- cbind(out, "kmer" = kmer.val[i])
    })
    out <- rbindlist(out)
    out[, lower.end := NULL]

    out <- as_tibble(out) %>% 
        mutate(kmer = paste0("kmer_", kmer)) %>% 
        tidyr::spread(key = "kmer", value = "upper.end")

    fwrite(
        x = out, 
        file = paste0("../data/Ranges_cutoffs.csv")
    )
}

# count number of curves to fit per exp
curve.counts <- df %>% 
    group_by(exp) %>%
    summarise(Count = n())

fwrite(
    x = curve.counts, 
    file = paste0("../data/kmer_", k, "_curve_counts.csv")
)

# Reassign org_file with curve counts per kmer
if(k == 8){
    files.to.import <- list.files(
        path = "../data/", 
        pattern = "curve_counts"
    )
    out <- lapply(1:length(files.to.import), function(i){
        out <- fread(paste0("../data/", files.to.import[i]))
        out <- cbind(out, "kmer" = kmer.val[i])
    })

    out <- rbindlist(out)

    out <- as_tibble(out) %>% 
        mutate(kmer = paste0("kmer_", kmer)) %>% 
        tidyr::spread(key = "kmer", value = "Count")

    org.file <- fread("../../Raw_data/org_file.csv")

    org.file <- as_tibble(org.file) %>% 
        mutate(
            exp = paste0(`Fragmentation type`, "/", `Experiment folder`), 
            .before = 1
        ) %>% 
        left_join(., out, by = "exp")

    org.file <- org.file %>% 
        dplyr::select(-exp)

    fwrite(
        x = org.file, 
        file = "../../Raw_data/org_file.csv"
    )
}