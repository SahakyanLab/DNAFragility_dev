# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.type <- as.character(args[2])
exp <- as.character(args[3])

# my.path = "/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/00_Preprocessing/scripts/"
# breakpoint.type="26-Chromosomal_fragile_sites"
# exp="5-Azacytidine"
setwd(my.path)

suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(library(data.table))

path.to.folder <- paste0("../../Raw_data/", breakpoint.type, "/", exp)
SOURCE <- stringr::str_extract(string = file, pattern = "^[^\\.]+")

df <- lapply(1:22, function(x){
    df <- plyranges::read_gff(paste0(path.to.folder, "/chr", x, "_fragile_gandharva.gff3"))
    df <- as_tibble(df)
    df <- df[, c(1,2,3,10,11)]
    colnames(df) <- c("seqnames", "start", "end", "name", "source")
    return(df)
})

df <- do.call(rbind, df)
df <- df %>% 
    mutate(
        across(where(is.factor), as.character),
        source = stringr::str_to_title(source)
    ) %>% 
    filter(source == exp) %>% 
    select(-source)

fwrite(
    x = df,
    file = paste0(path.to.folder, "/", exp, "_cfs_with_names.csv")
)

# save CFS region for kmertone
df <- df %>% 
    select(-name) %>% 
    rename_with(~c("chromosome","start.pos", "end")) %>% 
    mutate(strand = "+")

dir.create(paste0(path.to.folder, "/kmertone"), showWarnings = FALSE)

for(i in 1:22){
    cat("Extracting positions for chr", i, "...\n")

    filtered.df <- df %>% 
        filter(chromosome == paste0("chr", i))

    if(nrow(filtered.df) > 0){
        extracted.values <- lapply(1:nrow(filtered.df), function(x){
            seq(filtered.df$start.pos[x], filtered.df$end[x]-1, 1)
        })
        extracted.values <- unlist(extracted.values)
        extracted.values <- unique(extracted.values)

        DT <- data.table(
            chromosome = rep(paste0("chr", i), length(extracted.values)),
            start.pos = extracted.values,
            strand = "+"
        )

        fwrite(
            x = DT,
            file = paste0(path.to.folder, "/kmertone/chr", i, ".csv")
        )
    }
}