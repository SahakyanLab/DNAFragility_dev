suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(utils))

# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
reads_path <- as.character(args[1])
ref_path <- as.character(args[2])
my_path <- as.character(args[3])

setwd(paste0(my_path, "scripts/"))

# ---------------
pb <- txtProgressBar(min = 1, max = 23, style = 3)

# obtain breakpoint frequencies
for(i in 1:22){
  if(i == 23){
    chromosome = "chrX"
  } else if (i == 24) {
    chromosome = "chrY"
  } else {
    chromosome = paste0("chr",i)
  }
  
  # load data sets
  files <- list.files(path = paste0("../data/reads/", reads_path, "/breakpoint_positions/",
                                    chromosome, "/"),
                      pattern = "alignment_file_")

  tables <- lapply(files, function(x){
    read.table(paste0("../data/reads/", reads_path, "/breakpoint_positions/",
                      chromosome,"/", x),
               sep = ",", header = TRUE)
  })

  df <- do.call(rbind , tables) %>%
    as_tibble()
  colnames(df) <- c("start.pos", "lev.dist", "freq")

  # convert any non-integer columns to integers
  df <- df %>%
    mutate(across(where(is.character), as.integer)) %>%
    suppressWarnings()

  # remove any NAs
  df <- df %>%
    drop_na()

  # remove zero indices
  if(any(df$start.pos == 0)){
    df <- df %>% slice(which(df$start.pos != 0))
  }

  # save bottom ~95% of the levenshtein distance score
  # import average lev dist
  lev.dist.df <- read.table(file = paste0("../data/two_mers/", reads_path, "/AvgLevenshteinDistance.csv"),
                            sep = ",",
                            header = TRUE) %>% as_tibble()

  df <- df %>%
    filter(lev.dist < mean(lev.dist.df$SD1))

  # load reference sequence
  ref.seq <- readBStringSet(filepath = paste0("../data/ref/", ref_path, "/",
                                              chromosome, ".fasta"),
                            format = "fasta")
  # load dyad from reference sequence
  ref.two.mer <- read.csv(file = paste0("../data/two_mers/", reads_path, "/",
                                        chromosome, "_ref.csv"),
                          header = TRUE)

  # obtain all possible (4**k) k-mer permutations
  bases <- c("A", "G", "C", "T")
  two.mers <- permutations(n = length(bases),
                          v = bases,
                          r = 2,
                          repeats.allowed = TRUE)
  two.mers <- sapply(1:dim(two.mers)[1], function(x){
    paste(two.mers[x, ], collapse = "")
  })

  # account loaded dyad for strand symmetry
  rev.comp <- as.character(reverseComplement(DNAStringSet(two.mers)))
  two.mer.ref <- data.frame('fwd' = two.mers,
                            'rev.comp' = rev.comp) %>%
    mutate(cond = ifelse(seq(1:nrow(.)) < match(fwd, rev.comp), TRUE,
                         ifelse(fwd == rev.comp, TRUE, FALSE))) %>%
    filter(cond == TRUE) %>%
    select(-cond)

  breakpoints <- function(ind = 0){
    if(ind > 0){
      df <- df %>%
        mutate(start.pos = start.pos+ind)
    }

    # obtain indices
    fwd.ind <- match(two.mer.ref$fwd, ref.two.mer$value)
    rev.comp.ind <- match(two.mer.ref$rev.comp, ref.two.mer$value)

    # dyad count of the reference sequence
    two.mer.ref <- two.mer.ref %>%
      mutate(count = ref.two.mer$count[fwd.ind]+ref.two.mer$count[rev.comp.ind])

    # obtain dyads from alignment data
    mat <- matrix(data = c(df$start.pos,
                           df$start.pos+1),
                  ncol = 2)
    mat.ir <- IRanges(start = mat[,1], end = mat[,2])

    two.mer.seq <- paste(unlist(extractAt(x = ref.seq, at = mat.ir))) %>%
      as_tibble() %>%
      mutate(freq = df$freq) %>%
      group_by(value) %>%
      summarise(count = sum(freq)) %>%
      arrange(value) %>%
      filter(!str_detect(string = value, pattern = "N")) %>%
      rename(fwd = value)

    # final dyad frequency count
    two.mer.seq <- two.mer.ref %>%
      mutate(count = 0,
             count = two.mer.seq$count[fwd.ind]+two.mer.seq$count[rev.comp.ind])

    # breakage probability normalisation count
    breakpoints.df <- two.mer.seq %>%
      mutate(freq = count/two.mer.ref$count) %>%
      select(-count) %>%
      mutate(freq = freq/sum(freq))

    # return list
    return(list(breakpoints.df, two.mer.seq))
  }

  # save breakpoint data
  bp.probs <- breakpoints(ind = 0)

  bp.probs[[1]] %>%
    write.table(paste0("../data/two_mers/", reads_path, "/", chromosome, "_probs.txt"),
                sep = ",", row.names = FALSE)

  bp.probs[[2]] %>%
    write.table(paste0("../data/two_mers/", reads_path, "/", chromosome, "_freq.txt"),
                sep = ",", row.names = FALSE)

  setTxtProgressBar(pb, i)
}

# concatenate results
concat.breakpoints <- function(filepath){
  # load data sets
  files <- list.files(path = filepath,
                      pattern = "^chr([0-9]|1[0-9]|2[0-2]).*\\_freq.txt$")

  tables <- lapply(files, function(x){
    read.table(paste0(filepath, x),
               sep = ",", header = TRUE)
  })

  # sum of all frequencies
  freq.table <- do.call(cbind, lapply(tables, function(x) {
    rowSums(x[sapply(x, is.numeric)], na.rm = TRUE)
  }))
  freq.table <- rowSums(freq.table)

  # obtain cumulative kmer breakpoint frequencies
  kmer.freq <- tables[[1]] %>%
    as_tibble() %>%
    mutate(count = freq.table)

  # load dyad from reference sequence
  files <- list.files(path = filepath,
             pattern = "^chr([0-9]|1[0-9]|2[0-2]).*\\_ref.csv")

  tables <- lapply(files, function(x){
    read.table(paste0(filepath, x),
               sep = ",", header = TRUE)
  })

  freq.table <- do.call(cbind, lapply(tables, function(x){
    rowSums(x[sapply(x, is.numeric)], na.rm = TRUE)
  }))
  freq.table <- rowSums(freq.table)

  # obtain cumulative reference kmer breakpoint frequencies
  ref.freq <- tables[[1]] %>%
    as_tibble() %>%
    mutate(count = freq.table)

  # obtain all possible (4**k) k-mer permutations
  bases <- c("A", "G", "C", "T")
  two.mers <- permutations(n = length(bases),
                           v = bases,
                           r = 2,
                           repeats.allowed = TRUE)
  two.mers <- sapply(1:dim(two.mers)[1], function(x){
    paste(two.mers[x, ], collapse = "")
  })

  # account loaded dyas for strand symmetry
  rev.comp <- as.character(reverseComplement(DNAStringSet(two.mers)))
  two.mer.ref <- data.frame('fwd' = two.mers,
                            'rev.comp' = rev.comp) %>%
    mutate(cond = ifelse(seq(1:nrow(.)) < match(fwd, rev.comp), TRUE,
                         ifelse(fwd == rev.comp, TRUE, FALSE))) %>%
    filter(cond == TRUE) %>%
    select(-cond)

  # obtain indices
  fwd.ind <- match(two.mer.ref$fwd, ref.freq$value)
  rev.comp.ind <- match(two.mer.ref$rev.comp, ref.freq$value)

  # dyad count of the reference sequence
  ref.freq <- two.mer.ref %>%
    mutate(count = ref.freq$count[fwd.ind]+ref.freq$count[rev.comp.ind]) %>%
    as_tibble()

  # obtain final normalised breakpoint frequencies
  df <- kmer.freq %>%
    mutate(count = count/ref.freq$count) %>%
    mutate(count = count/sum(count))

  write.table(x = df, file = paste0("../data/two_mers/", reads_path, "/NormalisedBreakpointProb.csv"),
              sep = ",", row.names = FALSE)
  return(df)
}

df <- concat.breakpoints(filepath = paste0("../data/two_mers/", reads_path, "/"))
