setwd("/Volumes/Paddy_5TB/ProjectBoard_Patrick/03-Raw_Reads_Analysis/01_RMSD/scripts/")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))

#-----------------------------

RMSD <- function(a, b){
  round(sqrt(sum((a-b)^2)/length(a)))
}

#-----------------------------
# load data sets
chromosome = "chr1"

files <- list.files(path = "../../data/rmsd/exp_1/", 
                    pattern = "new_")
file.name <- paste0("kmer_", seq(from = 2, to = 8, by = 2))

# sort files
files <- str_sort(files, numeric = TRUE)

data.sets <- sapply(files, function(x){
  readRDS(file = paste0("../../data/rmsd/exp_1/", x))
})

# normalise by z-score where AUC = 1
# data.sets <- apply(data.sets, 2, scale, center = FALSE)
# colnames(data.sets) <- file.name

rmsd.plot <- data.sets %>%
  as_tibble() %>%
  mutate(id = -300.5:(nrow(data.sets)-300.5-1)) %>%
  melt(., id = "id") %>%
  ggplot(aes(x = id,
             y = value,
             col = variable)) + 
  geom_line(size = 0.8) +
  labs(x = "Position away from Breakpoint",
       y = "RMSD",
       title = "RMSD kmer frequencies")
rmsd.plot

ggsave(filename = paste0("../../figures/Simons_exp_1/", chromosome, "nofreq_RMSD.png"),
       plot = rmsd.plot)

plots <- lapply(1:dim(data.sets)[2], function(x){
  plots <- data.sets[,x] %>%
    as_tibble() %>%
    mutate(id = -300.5:(nrow(data.sets)-300.5-1)) %>%
    melt(., id = "id") %>%
    ggplot(aes(x = id,
               y = value)) + 
    geom_text(data = data.frame(
      xpos = Inf, 
      ypos =  Inf, 
      annotateText = file.name[x],
      hjustvar = 1 , 
      vjustvar = 1),
      aes(x = xpos,
          y = ypos,
          hjust = hjustvar,
          vjust = vjustvar,
          label = annotateText),
      size = 5) + 
    geom_line(size = 0.8) +
    labs(x = ifelse(x == 3, "Position away from Breakpoint", ""),
         y = ifelse(x == 3, "RMSD", ""))
  
  # return plots
  return(plots)
})

ggsave(filename = paste0("../../figures/Simons_exp_1/", chromosome, "sidebyside_nofreq_RMSD.png"),
       plot = do.call(grid.arrange, plots))
