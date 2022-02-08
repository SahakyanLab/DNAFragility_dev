# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
experiment <- as.character(args[2])
chromosome <- as.character(args[3])
setwd(paste0(my.path, "../lib/"))

my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/01_RMSD/scripts/"
experiment="03-Ancient_DNA/Altai_Neandertal_1"
# experiment="00-Ultrasonication/1000_Genomes_exp_1"
# experiment="01-Nebulization/1000_Genomes_exp_1"
# experiment="02-Sonication/1000_Genomes_exp_1"
chromosome="chr1"
setwd(paste0(my.path, "../lib/"))

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(stringr))

#-----------------------------
# load data sets
files <- list.files(path = paste0("../data/", experiment),
                    pattern = "rmsd_kmer_")
kmer.seq <- seq(from = 2, to = 10, by = 2)
file.name <- paste0("kmer_", kmer.seq)

# sort files
files <- str_sort(files, numeric = TRUE)

data.sets <- lapply(files, function(x){
  readRDS(file = paste0("../data/", experiment, "/", x))
})

P.value.plots <- lapply(1:length(file.name), function(x){
  rmsd.df <- data.sets[[x]]
  rmsd.df <- rmsd.df[301:length(rmsd.df)]
  Mean <- mean(rmsd.df)
  SD <- sd(rmsd.df)
  
  # noise computation
  noise <- rmsd.df[200:length(rmsd.df)]
  noise.mean <- mean(noise)
  noise.sd <- sd(noise)
  
  # calculate p-values
  p.value <- pnorm(rmsd.df,
                   mean = noise.mean,
                   sd = noise.sd,
                   lower.tail = T)
  tol <- -log10(0.05)
  sig.value <- -log10(p.value)
  
  # concatenate results
  df <- data.frame(rmsd = rmsd.df,
                   sig.value = sig.value) %>%
    as_tibble() %>%
    mutate(id = 0:(nrow(.)-1))
  
  # plot results
  p.value.plot <- df %>%
    ggplot(aes(x = id,
               y = sig.value)) + 
    annotate("rect", 
             xmin = 0, xmax = length(sig.value), 
             ymin = tol, ymax = ceiling(max(sig.value)),
             alpha = .7, fill = "lightgreen") +
    geom_line() + 
    ylim(0, ceiling(max(sig.value))) + 
    labs(x = "",
         y = ifelse(x == 1, "Sig Value",""),
         title = file.name[x]) + 
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # RMSD plot
  rmsd.plot <- df %>%
    ggplot(aes(x = id,
               y = rmsd)) + 
    geom_line() +
    geom_vline(xintercept = which(sig.value >= tol)[1],
               color = "darkgreen") + 
    geom_text(data = data.frame(xpos = which(sig.value > tol)[1]+1, 
                                ypos = max(df$rmsd),
                                annotateText = which(sig.value > tol)[1]-1,
                                hjustvar = 0, 
                                vjustvar = 1.1), 
              aes(x = xpos, 
                  y = ypos, 
                  hjust = hjustvar, 
                  vjust = vjustvar, 
                  label = annotateText),
              fontface = "bold", 
              size = 4) + 
    theme_bw() + 
    labs(x = ifelse(x == 1, "Position away from Breakpoint", ""),
         y = ifelse(x == 1, "RMSD", ""))
  
  # return plots
  g <- rbind(ggplotGrob(p.value.plot),
             ggplotGrob(rmsd.plot),
             size = "last")
  return(list(g, which(sig.value > tol)[1]+1))
})

g.plot <- grid.arrange(P.value.plots[[1]][[1]],
                       P.value.plots[[2]][[1]], 
                       P.value.plots[[3]][[1]], 
                       P.value.plots[[4]][[1]],
                       P.value.plots[[5]][[1]],
                       ncol = 5)

ggsave(filename = paste0("../figures/", experiment, "/", chromosome, "_RMSD_p-values.pdf"), 
       plot = g.plot, height = 7, width = 14)

# save threshold value for each k-mer
data.frame(
  kmer = kmer.seq,
  threshold = sapply(P.value.plots, "[[", 2)) %>%
  data.table::fwrite(file = paste0("../data/", experiment, "/", chromosome, "_threshold.csv"), 
                     row.names = FALSE)