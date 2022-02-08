# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
ref.seq <- as.character(args[3])
experiment.num <- as.character(args[4])
lower.limit <- as.integer(args[5])
upper.limit <- as.integer(args[6])
cores <- as.integer(args[7])

# my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/02_Correlations/scripts/"
# breakpoint.experiment="00-Ultrasonication/Simons_exp"
# experiment.num=1
# ref.seq="Simons_exp"
# lower.limit=1
# upper.limit=3
# cores=1
setwd(paste0(my.path, "../lib/"))

# load dependencies
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(pbapply))
source("../../00_ReadsAlignment/lib/LoadBreakpoints.R")
pbo = pboptions(type="txt")
# data.table::setDTthreads(threads = 1)

# iteration variables
chromosomes <- 1:22
exp <- seq(from = lower.limit, to = upper.limit, by = 1)
exp <- exp[-experiment.num]

# study one
cat("Loading in breakpoints from experiment", experiment.num, "...")
exp_1 <- lapply(
  chromosomes,
  LoadBreakpoints,
  path.to.origin = "../../", 
  experiment.folder = breakpoint.experiment, 
  experiment = experiment.num,
  select.col = "start.pos"
)
cat("Done!", "\n")

# overlaps between experiments
all.overlaps <- pbsapply(exp, function(e){
  overlaps.two.exp <- pbsapply(chromosomes, function(chr){
    # study 2
    exp_2 <- LoadBreakpoints(
      path.to.origin = "../../", 
      experiment.folder = breakpoint.experiment, 
      experiment = e,
      chromosome = chr,
      select.col = "start.pos"
    )
    
    # calculate overlap of breakpoints across experiments
    a <- exp_1[[chr]]
    b <- exp_2$start.pos
    
    overlap <- intersect(a, b)
    only.a <- a[-overlap]
    only.b <- b[-overlap]
    
    # absolute comparison
    overlap <- length(overlap)
    only.a <- length(only.a)
    only.b <- length(only.b)
    total <- overlap+only.a+only.b
    
    # percent comparison
    only.a.percent <- only.a/total
    only.b.percent <- only.b/total
    overlap.percent <- overlap/total
    return(overlap.percent)
  })
  return(overlaps.two.exp)
}, cl = cores)

colnames(all.overlaps) <- paste0("exp_", exp)

# all chromosomes
as.data.frame(all.overlaps) %>%
  fwrite(
    file = paste0("../data/overlap/", 
                  str_extract(string = breakpoint.experiment, pattern = "[0-9].+(?=/)"),
                  "_all_chromosomes.csv"),
    row.names = FALSE
  )