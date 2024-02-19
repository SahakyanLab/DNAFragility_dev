suppressPackageStartupMessages(suppressWarnings(library(data.table)))

download.file(
    "https://github.com/SahakyanLab/G4Damage/raw/master/raw_data/PQSdata.txt.gz", 
    "./PQSdata.txt.gz"
)
system("gunzip PQSdata.txt.gz")

df <- fread("./PQSdata.txt")
df <- df[, .(chr, str, genomic.start, genomic.end)]
setcolorder(df, c("chr", "genomic.start", "genomic.end", "str"))
setnames(df, c("Chromosome", "Start", "End", "Strand"))

# only keep autosomes
chr.to.keep <- paste0("chr", 1:22, "$", collapse = "|")
chr.to.keep.ind <- grepl(pattern = chr.to.keep, x = df$Chromosome)
chr.to.keep.ind <- which(chr.to.keep.ind)
df <- df[chr.to.keep.ind, ]

# sort custom
setorder(df, Start)
df <- df[order(match(Chromosome, paste0("chr", 1:22)))]

# save result
dir.create(
    path = "./kmertone",
    showWarnings = FALSE
)

# save chromosome separated files
for(chr in unique(df$Chromosome)){
    fwrite(
        df[Chromosome == chr],
        file = paste0("./kmertone/", chr, ".csv")
    )
}

fwrite(
    df,
    file = "./PQSdata.csv"
)

unlink("./PQSdata.txt")