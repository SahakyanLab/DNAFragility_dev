setwd("/Volumes/Paddy_5TB/ProjectBoard_Patrick/03-Raw_Reads_Analysis/scripts/ReadsAlignment/")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gtools))

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
  rename(start.pos = bp_start_pos,
         lev.dist  = lev_dist) %>%
  mutate(end.pos = start.pos+1) %>%
  relocate(end.pos, .after = start.pos) %>%
  mutate(across(everything(), as.numeric)) %>%
  drop_na()

group.by.neighbour <- function(tib){
  # Extract boolean values and cleans NA
  group.booleans <- paste(as.integer(tib$in.group), collapse = "")  
  group.booleans <- substr(group.booleans, 1, nchar(group.booleans)-1)
  # Split into groups at 0 delimiter
  groups <- unlist(str_split(group.booleans, "0"))
  # Apply function 
  groups.id <- sapply(1:length(groups), FUN = function(x){
    return(rep(x, nchar(groups[x])+1))
  })
  groups.id <- unlist(groups.id)
  return(groups.id)
}

normalise <- function(x){
  (x-min(x))/(max(x)-min(x))
}

frequency.mean <- function(tib){
  sum(tib$Diff*tib$freq)/sum(tib$freq) %>%
    suppressMessages()
}

frequency.sd <- function(tib, Mean){
  sqrt(sum((tib$Diff-Mean)**2*tib$freq)/(sum(tib$freq)-1)) %>%
    suppressMessages()
}

# save bottom ~95% of the levenshtein distance score
# import average lev dist
lev.dist.df <- read.table(file = "../../data/two_mers/Simons_exp_1/AvgLevenshteinDistance.csv",
                          sep = ",",header = TRUE) %>% as_tibble()

df <- df %>%
  filter(lev.dist < mean(lev.dist.df$SD2))

# double check for any NAs
which(is.na(df$start.pos))

# group adjacent breakpoint locations together
# group.id column will group together the breakpoints within gap.length
df.gaps <- function(gap.length = 1){
  gap.length <- seq(1, ceiling(gap.length))
  y <- df %>%
    mutate(Lag.Diff = (lag(start.pos)-start.pos)*-1,
           Next.Diff = lead(start.pos)-start.pos,
           in.group = Next.Diff %in% gap.length) %>%
    mutate(group.id = group.by.neighbour(.))
  return(y)
}

gap.length = 1
df.gaps.run <- df.gaps(gap.length = gap.length)

if(gap.length == 1){
  # frequency of adjacent breakpoint islands; gaps of length 1
  bp.freq <- df.gaps.run %>%
    group_by(group.id) %>%
    summarise(Diff = length(group.id)) %>%
    group_by(Diff) %>%
    summarise(freq = sum(Diff))
} else {
  # frequency of adjacent breakpoint islands; gaps of length > 1
  bp.freq <- df.gaps.run %>%
    group_by(group.id) %>%
    summarise(Diff.prelim = sum(diff(start.pos))) %>%
    mutate(Diff = ifelse(Diff.prelim == 0, 1, Diff.prelim)) %>%
    select(-Diff.prelim) %>%
    group_by(Diff) %>%
    summarise(freq = n()) 
}

# frequency distribution of length of consecutive breakpoints
bp.freq
bp.freq.plot <- bp.freq %>%
  filter(Diff <= 10) %>%
  ggplot(aes(Diff, freq)) +
  geom_bar(stat = "identity", fill = "#364f6b") +
  geom_text(aes(label = freq),
            position = position_dodge(width = 0.9), vjust = -0.3) +
  scale_x_continuous(breaks = 1:10,
                     labels = as.character(1:10)) +
  labs(x = "Length of Consecutive Breakpoints",
       y = "Frequency",
       title = paste0("Consecutive Breakpoints with Gap Length of: ", gap.length,"\n",
                      "Mean:    ", round(frequency.mean(bp.freq), 
                                         digits = 3),"\n",
                      "St. Dev: ", round(frequency.sd(bp.freq, 
                                                      frequency.mean(bp.freq)),
                                         digits = 3))) 

ggsave(paste0("../../figures/RawReadsIllumina/chr1/ReadDistribution_ConsecLength_",
              gap.length, ".pdf"),
       plot = bp.freq.plot,
       width = 12, 
       height = 8)

# frequency of the length of the gaps in-between apparent breakpoint islands
gaps.between.bp.freq <- df.gaps.run %>%
  group_by(group.id) %>%
  filter(row_number() %in% c(1, n())) %>%
  select(start.pos, group.id) %>%
  mutate(end.pos = start.pos[row_number() == n()]) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(Diff = lead(start.pos)-end.pos)

gaps.between.bp.freq <- gaps.between.bp.freq %>%
  group_by(Diff) %>%
  summarise(freq = n()) %>%
  na.omit() %>%
  filter(Diff > 0) %>%
  arrange(desc(freq)) 

# length of gaps between apparent breakpoint islands
gaps.between.bp.freq.plot <- gaps.between.bp.freq %>%
  top_n(n = 30) %>%
  ggplot(aes(Diff, freq)) +
  geom_bar(stat = "identity", fill = "#364f6b") +
  labs(x = "Length of gap between apparent breakpoint islands",
       y = "Frequency",
       title = paste0("Length of gap between apparent breakpoint islands\n",
                      "Mean:    ", round(frequency.mean(gaps.between.bp.freq), 
                                         digits = 3),"\n",
                      "St. Dev: ", round(frequency.sd(gaps.between.bp.freq, 
                                                      frequency.mean(gaps.between.bp.freq)),
                                         digits = 3)))

ggsave(paste0("../figures/RawReadsIllumina/chr1/ReadDistribution_GapBetweenBreakpointIslands.pdf"),
       plot = gaps.between.bp.freq.plot,
       width = 12, 
       height = 8)
