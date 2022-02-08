suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(pbapply))
pbo = pboptions(type="txt")
# data.table::setDTthreads(threads = 1)

# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
experiment.num <- as.character(args[3])
cores <- as.numeric(args[4])
setwd(paste0(my.path, "../lib"))

ConcatBreakpoints <- pblapply(1:22, function(i){
    # load data sets
    files <- list.files(
        path = paste0("../data/reads/", 
                      breakpoint.experiment, "_", experiment.num,
                      "/breakpoint_positions/chr", i, "/"),
        pattern = "alignment_file")
    
    files <- str_sort(files, numeric = TRUE)
    
    # concatenate files
    df <- lapply(files, function(x){
        fread(
            file = paste0("../data/reads/", 
                          breakpoint.experiment, "_", experiment.num,
                          "/breakpoint_positions/chr",
                          i, "/", x), 
            sep = ",",
            header = TRUE, 
            showProgress = FALSE)
    })

    df <- rbindlist(df)
    setnames(df, c("start.pos", "lev.dist", "freq"))
    
    # save bottom ~95% of the levenshtein distance score
    lev.dist.df <- fread(
        file = paste0("../data/average_levdist/", 
                      breakpoint.experiment, "_", experiment.num,
                      "/AvgLevenshteinDistance.csv"),
        sep = ",", 
        header = TRUE)

    df <- df[lev.dist < mean(lev.dist.df$SD1)]

    fwrite(
        df,
        file = paste0("../data/reads/", 
                      breakpoint.experiment, "_", experiment.num,
                      "/breakpoint_positions/chr", i, ".csv"),
        row.names = FALSE
    )
}, cl = cores)