# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
cores <- as.numeric(args[3])
max.chr <- as.numeric(args[4])

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(pbapply))
pbo = pboptions(type="txt")
# data.table::setDTthreads(threads = 1)

ConcatBreakpoints <- pblapply(1:max.chr, function(i){
    # load data sets
    files <- list.files(
        path = paste0("../../Raw_data/", breakpoint.experiment, 
                      "/breakpoint_positions/chr", i, "/"),
        pattern = "alignment_file")
    
    files <- str_sort(files, numeric = TRUE)
    
    # concatenate files
    df <- lapply(files, function(x){
        fread(
            file = paste0("../../Raw_data/", breakpoint.experiment,
                          "/breakpoint_positions/chr",
                          i, "/", x), 
            sep = ",",
            header = TRUE, 
            showProgress = FALSE)
    })

    df <- rbindlist(df)
    setnames(df, c("start.pos", "lev.dist", "freq"))
    
    # save bottom ~95% of the levenshtein distance score
    lev.dist.file <- paste0("../../Raw_data/", breakpoint.experiment, 
                            "/average_levdist/AvgLevenshteinDistance.csv")

    if(file.exists(lev.dist.file)){
        lev.dist.df <- fread(file = lev.dist.file, sep = ",", header = TRUE)
        df <- df[lev.dist < mean(lev.dist.df$SD1)]
    } else {
        df <- df[lev.dist < 1]
    }

    fwrite(
        df,
        file = paste0("../../Raw_data/", breakpoint.experiment, 
                      "/breakpoint_positions/chr", i, ".csv"),
        row.names = FALSE
    )
}, cl = cores)