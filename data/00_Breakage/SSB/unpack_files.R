suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))

# list folders
dirs <- c(
    "G5/6h", "G5/48", 
    "Imatinib/6h", "Imatinib/48h", 
    "SN/6h", "SN/48h", 
    "Talazoparib/6h", "Talazoparib/48h"
)

# for each folder, concatenate files
for(x in 1:length(dirs)){
    files <- list.files(
        path = paste0("./", dirs[x]),
        pattern = ".bed",
        full.names = TRUE,
    )
    df <- lapply(files, plyranges::read_bed)
    df <- do.call(plyranges::bind_ranges, df)
    
    df <- df %>% 
        unique() %>% 
        dplyr::select(-name, -score)

    seqlevels(df, pruning.mode = "coarse") <- paste0("chr", 1:22)
    df <- df %>% 
        dplyr::arrange(seqnames, start)

    plyranges::write_bed(
        df,
        paste0("./", dirs[x], "/output.bed")
    )   

    # rm old files
    unlink(files, force = TRUE)
}