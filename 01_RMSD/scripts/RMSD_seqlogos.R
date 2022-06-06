my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/01_RMSD/scripts/"
k=4
setwd(my.path)

suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(seqLogo))
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
      file = paste0("../data/", file, "/freq_rmsd_kmer_", k, ".Rdata")
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
      count = out
    ) %>% 
    mutate(category = category)

  return(df.kmer)
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

# categories <- unique(org.file$Category)
categories <- c("Mechanical", "Biological", "Cell-free", "Ancient")

for(c in categories){
  df.kmer <- pblapply(categories, function(category){
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
  df.kmer <- do.call(rbind, df.kmer)

  # extract top/bottom percent by count
  get.subset <- function(which.category, percent){
      temp <- df.kmer %>% 
          filter(category == which.category)

      N <- ceiling(percent/100*nrow(temp))
      return(
          rbind(
              temp %>% slice_max(order_by = count, n = N),
              temp %>% slice_min(order_by = count, n = N)
          )
      )
  }

  percent=10
  df <- get.subset(which.category = c, percent = percent)
  df <- df %>% 
      mutate(count = count/n())

  enriched.kmers <- rep(df$kmer, df$count)
  m <- consensusMatrix(enriched.kmers, as.prob = TRUE)
  p <- seqLogo::makePWM(m)
  pdf(paste0(c, "-", k, "-kmer_seqlogo_info_bits_enriched.pdf")) 
  seqLogo(p, ic.scale = TRUE)
  plot.save <- dev.off()

  depleted.kmers <- paste0(reverseComplement(DNAStringSet(enriched.kmers)))
  m <- consensusMatrix(depleted.kmers, as.prob = TRUE)
  p <- seqLogo::makePWM(m)
  pdf(paste0(c, "-", k, "-kmer_seqlogo_info_bits_depleted.pdf")) 
  seqLogo(p, ic.scale = TRUE)
  plot.save <- dev.off()
}