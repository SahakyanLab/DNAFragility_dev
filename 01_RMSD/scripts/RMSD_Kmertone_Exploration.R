my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/01_RMSD/scripts/"
k=4
setwd(my.path)

suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(pbapply))
if(length(args) > 0) pbo = pboptions(type="txt")

source("../lib/CalcKmerFreq.R")
source("../lib/CalcRMSD.R")
source("../../00_ReadsAlignment/lib/LoadBreakpoints.R")
source("../../02_Correlations/lib/GenerateKmerTable.R")
source("../../02_Correlations/lib/LoadKmertoneData.R")
# data.table::setDTthreads(threads = 1)

# create reference table of kmers
sample.kmer <- data.table(kmer = GenerateKmerTable(k, ref = FALSE))
kmer.ref <- GenerateKmerTable(k, ref = TRUE)

# import range of interaction limits
df.ranges <- fread("../data/Ranges_cutoffs.csv")
hash_map <- as_tibble(df.ranges) %>% 
    select(Cluster, paste0("kmer_", k)) %>% 
    pull(paste0("kmer_", k), Cluster)

# filter for ranges
which.range <- "short.range"
each.dir <- floor(unname(hash_map[which.range]/2))
range.header <- str_replace(string = which.range, pattern = "\\.", replacement = "-")

get.data <- function(org, sample.kmer, kmer.ref, which.range, category){
  org.filtered <- org %>% 
    filter(Category == category) %>% 
    pull(breakpoint.experiment)

  # import RMSD results
  norm.kmer.count <- lapply(org.filtered, function(file){
    out <- readRDS(
      file = paste0("../data/", file, "/normalised_freq_rmsd_kmer_", k, ".Rdata")
    )

    # extract range of interest
    limits <- ncol(out)/2-1
    limits <- -limits:(ncol(out)-limits-1)

    lower.end <- which(limits %in% -which.range)
    upper.end <- which(limits %in% which.range)+1

    return(rowSums(out[, lower.end:upper.end]))
  })
  out <- do.call(cbind, norm.kmer.count)
  out <- rowMeans(out)
  df.kmer <- tibble(
      kmer = kmer.ref$kmer,
      norm.mean.count = out
    )

  # import kmertone results
  import.kmertone <- function(action, file){
    out <- LoadKmertoneData(
      breakpoint.experiment = file,
      k = k,
      kmer.table = sample.kmer$kmer,
      kmer.ref.table = kmer.ref,
      action = action,
      to.target.dir = "../../02_Correlations/"
    )
    return(out)
  }

  actions <- c("z-score", "ratio")
  out <- pbsapply(actions, function(action){
    out <- lapply(org.filtered, import.kmertone, action = action)

    out <- sapply(1:length(out[[1]]), function(x){
      mean(sapply(out, `[[`, x))
    })

    return(out)
  })

  df.kmertone <- as_tibble(out) %>% 
    mutate(
      kmer = kmer.ref$kmer, 
      .before = 1
    ) %>% 
    rename_with(~c("kmer", "z.score", "ratio"))

  df.rmsd.kmertone <- left_join(
    df.kmer, df.kmertone,
    by = "kmer"
  )

  df.rmsd.kmertone <- df.rmsd.kmertone %>% 
    mutate(category = category)

  return(df.rmsd.kmertone)
}

# import breakpoint data of specified group of breakpoint
org.file <- fread("../../Raw_data/org_file.csv")
org.file <- as_tibble(org.file) %>% 
  mutate(
    fragmentation.type = `Fragmentation type`,
    breakpoint.experiment = paste0(fragmentation.type, "/", `Experiment folder`),
    .before = 1
  ) %>% 
  select(-`Fragmentation type`)

dirs <- list.dirs(
  path = paste0("../data"),
  full.names = FALSE,
  recursive = FALSE
)

to.select <- !is.na(match(org.file$fragmentation.type, dirs))
org.file <- org.file[to.select, ]

categories <- unique(org.file$Category)

########################################################################
# org.file
categories <- c("Mechanical", "Biological", "Cell-free", "Ancient")


# org = org.file
# sample.kmer = sample.kmer 
# kmer.ref = kmer.ref 
# which.range = each.dir 
# category = categories[1]

df.rmsd.kmertone <- pblapply(categories, function(category){
  return(
    get.data(
      org = org.file, 
      sample.kmer = sample.kmer, 
      kmer.ref = kmer.ref, 
      which.range = each.dir, 
      category = category
    )
  )
})
df.rmsd.kmertone <- do.call(rbind, df.rmsd.kmertone)

# Identify which kmers are enriched and depleted
# within the specific short/mid/long-range of interaction
suppressPackageStartupMessages(library(ggplot2))

CI <- t.test(df.rmsd.kmertone$norm.mean.count)$conf.int

get.subset <- function(which.category, N){
  temp <- df.rmsd.kmertone %>% 
    filter(category == which.category)
  return(
    rbind(
      temp %>% slice_max(order_by = norm.mean.count, n = N),
      temp %>% slice_min(order_by = norm.mean.count, n = N)
    )
  )
}

N=5
df.subset <- lapply(categories, function(x){
  get.subset(which.category = x, N=N)
})
df.subset <- do.call(rbind, df.subset)

edge.val <- max(abs(log2(df.rmsd.kmertone$ratio)))
edge.val <- edge.val*1.02

p <- df.rmsd.kmertone %>% 
  ggplot(aes(
    x = log2(ratio),
    y = abs(z.score),
    # size = norm.mean.count,
    color = norm.mean.count
  )) + 
  geom_point(alpha = 0.4) +
  facet_grid(~category) + 
  scale_color_gradient(
    limits = CI,
    low = "#00318b", 
    high = "#ff3c00",
    oob = scales::squish
  ) + 
  ggrepel::geom_text_repel(
    data = df.subset, 
    aes(
      x = log2(ratio),
      y = abs(z.score),
      label = kmer
    ), 
    box.padding = 0.5, 
    max.overlaps = Inf,
    min.segment.length = 0,
    size = 3
  ) + 
  labs(
    title = paste0(k, "-mer contributing to ", range.header, " sequence influence"),
    x = "Log2 average probability ratios", 
    y = "Average absolute z-scores", 
    colour = "Average norm. \nk-mer freq"
  )

ggsave(
  plot = p,
  filename = paste0("../figures/", k, "-mers_contributing_", range.header, ".pdf"),
  height = 7, width = 13
)