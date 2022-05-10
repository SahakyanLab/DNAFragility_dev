# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
chromosome <- as.character(args[3])
category <- as.character(args[4])
setwd(my.path)

suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))

source("../lib/FitGMM.R")
source("../lib/IntegralFunctions.R")
source("../lib/CalcCI.R")
source("../lib/MakePlot.R")

# load data sets
files <- list.files(
  path = paste0("../data/", breakpoint.experiment),
  pattern = "rmsd_kmer_"
)
file.name <- paste0("kmer_", seq(from = 2, to = 10, by = 2))
files <- str_sort(files, numeric = TRUE)

data.sets <- lapply(1:length(files), function(x){
  out <- readRDS(
    file = paste0("../data/", breakpoint.experiment,"/", files[x])
  )

  limits <- length(out)/2-1
  out <- as_tibble(out) %>% 
    mutate(
      kmer = as.factor(file.name[x]), 
      x = -limits:(length(out)-limits-1)
    ) %>% 
    dplyr::rename(y = value)
  
  return(out)
})

data.sets <- do.call(rbind, data.sets)

# Overall RMSD plots
plots <- data.sets %>%
  ggplot(aes(x = x, y = y)) + 
  geom_line(size = 0.8) + 
  facet_wrap(~kmer, ncol = 2, scales = "free_y") + 
  theme_bw() + 
  labs(
    x = "Position away from breakpoint",
    y = "RMSD"
  )

ggsave(
  filename = paste0("../figures/", breakpoint.experiment, 
                    "/chr", chromosome, "_sidebyside_RMSD.pdf"),
  plot = plots, 
  height = 7, width = 7
)

# GMM plots
data.sets <- data.sets %>%
  mutate(across(where(is.factor), as.character)) %>% 
  filter(kmer != "kmer_2" & kmer != "kmer_10")

p <- vector(mode = "list", length = 3)

fits <- c("kmer_4" = 3, "kmer_6" = 3, "kmer_8" = 3)

for(k in 1:length(fits)){
  print(k)
  dat <- data.sets %>% filter(kmer == names(fits[k]))
  curvefits <- FitGMM(dat = dat, ind = names(fits[k]), nr.of.curves = fits[k])

  if(class(curvefits[length(curvefits)][[1]]) != "nls"){
    fits[k] <- 2
    curvefits <- FitGMM(dat = dat, ind = names(fits[k]), nr.of.curves = fits[k])

    if(class(curvefits[length(curvefits)][[1]]) != "nls"){
      fits[k] <- 1
      curvefits <- FitGMM(dat = dat, ind = names(fits[k]), nr.of.curves = fits[k])
    }

    if(all(is.na(unlist(curvefits[length(curvefits)-2])))){
      p <- p[1:(length(p)-1)]
      break
    }
  }

  output <- MakePlot(dat = dat, k = names(fits[k]), 
                     curve.vals = curvefits, 
                     nr.of.curves = fits[k])
  output.plot <- output[[1]]

  if(class(output.plot)[1] == "numeric"){
    fits[k] <- 2
    curvefits <- FitGMM(dat = dat, ind = names(fits[k]), nr.of.curves = fits[k])

    if(class(curvefits[length(curvefits)][[1]]) != "nls" & names(fits[k]) != "kmer_8"){
      fits[k] <- 1
      curvefits <- FitGMM(dat = dat, ind = names(fits[k]), nr.of.curves = fits[k])
    }

    if(all(is.na(unlist(curvefits[length(curvefits)-2])))){
      p <- p[1:(length(p)-1)]
      break
    }

    output <- MakePlot(dat = dat, k = names(fits[k]), 
                       curve.vals = curvefits, 
                       nr.of.curves = fits[k])
    output.plot <- output[[1]]
  } else {
    p[[k]] <- output.plot
  }

  # percent contribution of each gaussian curve towards breakage
  df <- as_tibble(curvefits[length(curvefits)-1][[1]])
  df[1,] <- df[1,]*100

  # range of influence based on 95 percent confidence intervals
  output.CIlst <- output[[2]]
  curves <- output.CIlst[seq(2, length(output.CIlst), 2)]-
            output.CIlst[seq(1, length(output.CIlst), 2)]

  # sd of curves
  SD <- sapply(1:(length(curvefits)-2), function(x){
    vals <- curvefits[[x]]

    left <- t.test(vals, conf.level = 0.68)$conf.int[1]
    right <- t.test(vals, conf.level = 0.68)$conf.int[2]

    lower.bound <- which(vals > left)[1]
    upper.bound <- which(vals > right)
    upper.bound <- upper.bound[length(upper.bound)]

    return((upper.bound-lower.bound)/2)
  })

  SD <- c(SD, rep(NA_real_, 3-length(SD)))

  # peak intensity of each gaussian curve
  peak.intensity <- sapply(1:fits[k], function(x){
    max(curvefits[[x]])
  })

  peak.intensity <- c(peak.intensity, rep(NA_real_, 3-length(peak.intensity)))

  # combine results
  df <- rbind(df, curves, SD, peak.intensity)
  df <- df %>% 
    mutate(
      exp = breakpoint.experiment,
      category = category,
      rowid = c("contribution", "ranges", "SD", "peak.intensty"), 
      .before = 1,
    )

  write.csv(
    x = df,
    file = paste0("../data/", breakpoint.experiment, "/key_stats_", names(fits[k]), ".csv"),
    row.names = FALSE
  )

  # if(k == 1) print(p) # temp. avoid an unsolved ggplot bug
}

ggsave(
  filename = paste0("../figures/", breakpoint.experiment, 
                    "/chr", chromosome, "_RMSD_all_GMMFit.pdf"),
  plot = cowplot::plot_grid(plotlist = p, axis = "b", ncol = 1), 
  height = 15, width = 7
) 