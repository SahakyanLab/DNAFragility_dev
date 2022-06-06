# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.type <- as.character(args[2])
exp <- as.character(args[3])
file.format <- as.character(args[4])
start.idx <- as.numeric(args[5])
end.idx <- as.numeric(args[6])
alignment_strands <- as.character(args[7])
fastq.processed <- as.logical(args[8])

# my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/00_Preprocessing/scripts/"

# breakpoint.type="10-DSBCapture"
# exp="NHEK_DSBs/archive/NHEK_DSBs"

# file.format="BED"
# start.idx=-1
# end.idx=0
# alignment_strands="both"

# fastq.processed=TRUE
setwd(my.path)

# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(stringr))
source("../lib/LoadFiles.R")
source("../lib/ProcessFiles.R")

path.to.save.bp <- paste0("../../Raw_data/", breakpoint.type, "/", exp)

all.files <- list.files(
    path = paste0("../../Raw_data/", breakpoint.type, "/", exp),
    full.names = TRUE,
)

files.to.process <- all.files[!file.info(all.files)$isdir]
files.to.process <- str_sort(files.to.process, numeric = TRUE)

do.not.process <- str_extract(
    string = files.to.process, 
    pattern = "fastq|fasta|bam|bai"
)
files.to.process <- files.to.process[is.na(do.not.process)]

file.names <- str_sort(str_extract(
    string = files.to.process, 
    pattern = "chr[[:alnum:]]+"), 
    numeric = TRUE
)

both.ends <- ifelse(alignment_strands == "both", TRUE, FALSE)

for(i in seq_along(files.to.process)){
    data <- LoadFiles(
        data = files.to.process[i], 
        file.type = file.format,
        fastq.processed = fastq.processed
    )

    if(all(grepl(pattern = "plus|minus", x = basename(files.to.process)))){
        if(grepl(pattern = "plus", x = basename(files.to.process)[i])){
            strands <- "plus"
        } else if(grepl(pattern = "minus", x = basename(files.to.process)[i])){
            strands <- "minus"
        }
    } else {
        strands <- ""
    }

    if(all(paste0("chr", 1:22) %in% file.names)){
        cat("Processing chr", i, "\n")
        ProcessFiles(
            data = data, 
            chr = str_extract(string = file.names[i], pattern = "[[:digit:]]+"),
            file.type = file.format,
            start.index = start.idx,
            end.index = end.idx,
            path.to.save.bp = path.to.save.bp,
            both.ends = both.ends,
            strands = strands,
            fastq.processed = fastq.processed
        )
    } else {
        for(chr in 1:22){
            cat("Processing chr", chr, "\n")
            ProcessFiles(
                data = data,
                chr = chr,
                file.type = file.format,
                start.index = start.idx,
                end.index = end.idx,
                path.to.save.bp = path.to.save.bp,
                both.ends = both.ends,
                strands = strands,
                fastq.processed = fastq.processed
            )
        }
    }
}

for(chr in 1:22){
    if(strands != ""){
        dt.plus <- fread(paste0(path.to.save.bp, "/breakpoint_positions/plus_chr", chr, ".csv"))
        dt.minus <- fread(paste0(path.to.save.bp, "/breakpoint_positions/minus_chr", chr, ".csv"))

        dt <- rbind(dt.plus, dt.minus)

        if(!"strand" %in% colnames(dt)){
            dt[, strand := "+"]    
        }
        
        setorder(dt, strand)
        fwrite(dt, paste0(path.to.save.bp, "/breakpoint_positions/chr", chr, ".csv"))

        system(paste0("rm ", path.to.save.bp, "/breakpoint_positions/plus_chr", chr, ".csv"))
        system(paste0("rm ", path.to.save.bp, "/breakpoint_positions/minus_chr", chr, ".csv"))
    }

    # save breakpoint positions for kmertone
    dt <- fread(paste0(path.to.save.bp, "/breakpoint_positions/chr", chr, ".csv"))
    dt[, chromosome := paste0("chr", chr)]

    if(!"strand" %in% colnames(dt)){
        dt[, strand := "+"]    
    }

    setcolorder(dt, c("chromosome", "start.pos", "strand"))
    fwrite(dt, paste0(path.to.save.bp, "/kmertone/chr", chr, ".csv"))
}