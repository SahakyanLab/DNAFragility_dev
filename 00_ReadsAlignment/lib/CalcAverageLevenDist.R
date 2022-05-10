# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
cores <- as.numeric(args[3])

# my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/00_ReadsAlignment/scripts/"
# breakpoint.experiment="00-Ultrasonication/Simons_exp_25"
# cores=1
setwd(paste0(my.path, "../lib"))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(stringr))
pbo = pboptions(type="txt")
# data.table::setDTthreads(threads = 1)

# ---------------
# obtain average levenshtein distance for all chromosomes
lev.dist <- pblapply(1:22, function(i){  
  files <- list.files(
    path = paste0("../../Raw_data/", breakpoint.experiment, 
                  "/breakpoint_positions/chr", i),
    pattern = "alignment_file"
  )
  files <- str_sort(files, numeric = TRUE)
  
  tables <- lapply(files, function(x){
    fread(
      paste0("../../Raw_data/", breakpoint.experiment, 
             "/breakpoint_positions/chr", i, "/", x), 
      sep = ",", 
      header = TRUE,
      drop = "bp.start.pos",
      showProgress = FALSE
    )
  })

  df <- rbindlist(tables)
  setnames(df, c("lev.dist", "freq"))

  if(anyNA(df)){
    df <- df[complete.cases(df)]
  }

  if("character" %in% sapply(df, class)){
    df[, names(df) := lapply(.SD, as.numeric)]
  }

  # return frequency average and st.dev levenshtein distance
  lev.dist.df <- df[, .(count = .N), by = lev.dist]
  Mean <- sum(lev.dist.df$lev.dist*lev.dist.df$count)/sum(lev.dist.df$count)
  SD <- sqrt(sum((lev.dist.df$lev.dist)**2*lev.dist.df$count)/(sum(lev.dist.df$count)-1))
  
  return(list(Mean, SD))
}, cl = cores)

# concatenate results and calculate overall statistics
df <- rbindlist(lev.dist)

if("list" %in% sapply(df, class)){
  df[, names(df) := lapply(.SD, as.numeric)]
}

setnames(df, c("Mean", "SD"))
df[, `:=`(SD1 = Mean+SD, SD2 = Mean+2*SD)]

# rename column
df[, Chromosomes := 1:nrow(df)]
setorder(df, -Chromosomes)
df[, Chromosomes := paste0("Chr", 22:1)]
df[, `:=`(Chromosomes = forcats::fct_inorder(Chromosomes))]

# obtain bar plot of levenshtein distance
lev.plot <- ggplot(df) +
  geom_bar(
    aes(x = Chromosomes,
        y = Mean),
    stat = "identity",
    fill = "skyblue",
    alpha = 0.7) + 
  geom_pointrange(
    aes(x = Chromosomes,
        y = Mean,
        ymin = Mean,
        ymax = SD1),
    color = "orange",
    alpha = 1,
    size = 0.8) +
  coord_flip() +
  labs(
    title = "Average Levenshtein Distance plotted as Mean + 1 St.Dev", 
    subtitle = paste0("Overall Average: ", signif(mean(df$Mean), 3)))

ggsave(filename = paste0("../figures/", breakpoint.experiment, 
                        "/AvgLevenshteinDistance.pdf"),
      plot = lev.plot)

# save lev dist 
fwrite(x = df, 
       file = paste0("../../Raw_data/", breakpoint.experiment, 
                     "/average_levdist/AvgLevenshteinDistance.csv"), 
       row.names = FALSE)