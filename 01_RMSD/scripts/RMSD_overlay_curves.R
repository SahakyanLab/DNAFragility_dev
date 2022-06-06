my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/01_RMSD/scripts/"
setwd(my.path)

suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pbapply))

source("../lib/CalcKmerFreq.R")
source("../lib/CalcRMSD.R")
source("../../00_ReadsAlignment/lib/LoadBreakpoints.R")
source("../../02_Correlations/lib/GenerateKmerTable.R")
source("../../02_Correlations/lib/LoadKmertoneData.R")
# data.table::setDTthreads(threads = 1)

# import breakpoint data of specified group of breakpoint
org.file <- fread("../../Raw_data/org_file.csv")
org.file <- as_tibble(org.file) %>% 
  mutate(
    fragmentation.type = `Fragmentation type`,
    breakpoint.experiment = paste0(fragmentation.type, "/", `Experiment folder`),
    .before = 1
  ) %>% 
  select(-`Fragmentation type`) %>% 
  filter(`DSB Map` == TRUE)

dirs <- list.dirs(
  path = paste0("../data"),
  full.names = FALSE,
  recursive = FALSE
)

to.select <- !is.na(match(org.file$fragmentation.type, dirs))
org.file <- org.file[to.select, ]
categories <- unique(org.file$Category)
categories <- c("Mechanical", "Biological", "Cell-free", "Ancient", "Enzymatic")

get.data <- function(org, sample.kmer, kmer.ref, category, k){
  org.filtered <- org %>% 
    filter(Category == category) %>% 
    pull(breakpoint.experiment)

  # import RMSD results
  import.rmsd <- lapply(org.filtered, function(file){
    out <- readRDS(
      file = paste0("../data/", file, "/rmsd_kmer_", k, ".Rdata")
    )
    
    limits <- length(out)/2-1
    out <- as_tibble(out) %>%
      mutate(
        x = -limits:(length(out)-limits-1),
        category = category
      ) %>% 
      dplyr::rename(y = value)

    return(out)
  })
  out <- do.call(rbind, import.rmsd)
  return(out)
}

k=8
out <- pblapply(categories, function(category){
  return(
    get.data(
      org = org.file, 
      sample.kmer = sample.kmer, 
      kmer.ref = kmer.ref, 
      category = category
    )
  )
})
df <- do.call(rbind, out)

average.curve <- df %>% 
  group_by(category, x) %>% 
  summarise(RMSD = mean(y))

p <- df %>% 
  ggplot(aes(
    x = x,
    y = y
  )) + 
  geom_line(alpha = 0.2) + 
  geom_line(
    data = average.curve,
    aes(
      x = x,
      y = RMSD
    )
  ) +
  facet_wrap(~category, scales = "free_y") + 
  labs(
    title = paste0("Kmer: ", k),
    x = "Position away from breakpoint",
    y = "RMSD"
  )

ggsave(
  plot = p,
  filename = paste0("../figures/", k, "-mers_overlaid-all-RMSD.pdf"),
  height = 7, width = 13
)

out <- pblapply(seq(2,8,2), function(k){
  out <- pblapply(categories, function(category){
    return(
      get.data(
        org = org.file, 
        sample.kmer = sample.kmer, 
        kmer.ref = kmer.ref, 
        category = category,
        k = k
      )
    )
  })
  df <- do.call(rbind, out)
  df <- df %>% 
    mutate(kmer = paste0("kmer_", k))
  return(df)
})
df <- do.call(rbind, out)

average.curve <- df %>% 
  group_by(category, x, kmer) %>% 
  summarise(RMSD = mean(y))

p <- average.curve %>% 
  ggplot(aes(
    x = x,
    y = RMSD
  )) + 
  geom_line() + 
  facet_grid(
    vars(kmer), 
    vars(category),
    scales = "free_y"
  ) + 
  labs(
    x = "Position away from breakpoint",
    y = "RMSD"
  )

ggsave(
  plot = p,
  filename = "../figures/Average-RMSD.pdf",
  height = 7, width = 10
)