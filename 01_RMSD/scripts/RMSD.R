setwd("/Volumes/Paddy_5TB/ProjectBoard_Patrick/03-Raw_Reads_Analysis/01_RMSD/scripts/")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(gtools))

#-----------------------------
# helper functions
RMSD <- function(a, b){
  return(sqrt(mean((a-b)^2, na.rm = TRUE)))
}

#-----------------------------
# load data sets
chromosome = "chr1"

files <- list.files(path = paste0("../../data/reads/Simons_exp_1/breakpoint_positions/",chromosome,"/"))

tables <- lapply(files, function(x){
  fread(paste0("../../data/reads/Simons_exp_1/breakpoint_positions/",chromosome,"/", x), 
        sep = ",", header = TRUE)
})

df <- do.call(rbind , tables) %>% 
  as_tibble() %>%
  rename_with(~c("start.pos", "lev.dist", "freq"))

# save bottom ~95% of the levenshtein distance score
# import average lev dist
lev.dist.df <- read.table(file = "../../data/two_mers/Simons_exp_1/AvgLevenshteinDistance.csv",
                          sep = ",", header = TRUE) %>% as_tibble()

df <- df %>%
  filter(lev.dist < mean(lev.dist.df$SD1)) %>%
  select(start.pos)

# double check for any NAs
which(is.na(df))

# load reference sequence
ref.seq <- readBStringSet(filepath = paste0("../../data/ref/Simons_exp/", chromosome, ".fasta"), 
                          format = "fasta")

iters <- seq(from = 2, to = 10, by = 2)
rmsd.all.values <- lapply(iters, function(k){
  # obtain all possible (4**k) k-mer permutations
  bases   <- c("A", "G", "C", "T")
  
  # generate all k-mers 
  k.mers <- permutations(n = length(bases), v = bases, 
                         r = k, repeats.allowed = TRUE)
  k.mers <- sapply(1:dim(k.mers)[1], function(i){
    paste(k.mers[i,], collapse = "")
  })
  
  # account loaded dyas for strand symmetry 
  rev.comp <- as.character(reverseComplement(DNAStringSet(k.mers)))
  k.mer.ref <- data.frame('fwd'      = k.mers,
                          'rev.comp' = rev.comp) %>%
    mutate(cond = ifelse(seq(1:nrow(.)) < match(fwd, rev.comp), 
                         TRUE, ifelse(fwd == rev.comp, TRUE, FALSE))) %>%
    filter(cond == TRUE) %>%
    select(-cond)
  
  freq.calc.i <- function(ind, k){
    # dyad frequency of iteration, i
    df <- df %>%
      mutate(start.pos = start.pos+ind)
    
    if ((k %% 2) == 0) {
      # even numbers
      end.pos = ceiling((k-1)/2)
      start.pos = (k-1)-end.pos
    }

    # obtain k-mers from alignment data
    mat <- matrix(data = c(df$start.pos-start.pos,
                           df$start.pos+end.pos), ncol = 2)
    
    # extract k-meric counts and relative frequencies
    k.mer.seq <- str_sub(string = ref.seq, start = mat[,1], end = mat[,2]) %>%
      as_tibble() %>%
      group_by(value) %>%
      summarise(n = n()) %>%
      filter(!str_detect(string = value, pattern = "N")) %>% 
      dplyr::rename(fwd = value)
    
    # account for strand symmetry
    fwd.ind <- match(k.mer.ref$fwd, k.mer.seq$fwd)
    rev.comp.ind <- match(k.mer.ref$rev.comp, k.mer.seq$fwd)
    
    # update data frame with dyad frequency count
    k.mer.i <- k.mer.ref %>%
      mutate(freq = k.mer.seq$n[fwd.ind]+k.mer.seq$n[rev.comp.ind],
             freq = freq/sum(freq, na.rm = TRUE))
    
    # return data frame of k-mers with associated breakpoint frequencies
    return(k.mer.i)
  }
  
  freq <- pbsapply(-300:301, function(x){
    freq.calc.i(ind = x, k = k)$freq
  })

  rmsd.values <- pbsapply(1:(dim(freq)[2]-1), function(x){
    RMSD(freq[,x], freq[,(x+1)])
  })
  
  # save output
  saveRDS(object = rmsd.values,
          file = paste0("../../data/rmsd/exp_1/new_rmsd_kmer_", k, ".Rdata"))
  
  # return(rmsd.values)
})
