# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)
the.category="Ultrasonication"

for(k in c(6,8)){
    # probability ratio
    df.ratio <- fread(paste0("../data/kmertone/QueryTable/", 
                              "QueryTable_kmer-", k, "_ratio.csv"))
    df.ratio <- df.ratio[category == the.category, -"category"]
    df.ratio <- t(df.ratio)
    df <- as.data.table(df.ratio)
    df[, kmer := rownames(df.ratio)]
    setnames(df, c("prob", "kmer"))
    setcolorder(df, c("kmer", "prob"))

    # extend table to have all k-mers, not just lexicographically in the correct order
    rev.comp <- as.character(
        Biostrings::reverseComplement(
            Biostrings::DNAStringSet(df$kmer)
        )
    )
    df[, rev_comp := rev.comp]

    k.mers <- do.call(
        data.table::CJ,
        rep(list(c("A", "C", "G", "T")), k)
    )
    kmer_list <- k.mers[, do.call(paste0, .SD)]
    kmer_ref <- data.table(kmer=kmer_list)
    kmer_ref[, prob := 0]

    # forward k-mers
    ind_fwd <- match(df$kmer, kmer_ref$kmer)    
    kmer_ref$prob[ind_fwd] <- df$prob

    # reverse complement k-mers
    ind_rc <- match(df$rev_comp, kmer_ref$kmer)
    ind_rc[is.na(ind_rc)] <- (nrow(kmer_ref)-nrow(df)+1):(nrow(kmer_ref))
    kmer_ref$prob[ind_rc] <- df$prob
    df[, rev_comp := NULL]

    # save results
    fwrite(
        kmer_ref, 
        file = paste0("../../../github_repos/GenomeAssembler_dev/data/", 
                      "QueryTable/QueryTable_kmer-", k,".csv"),
        showProgress = FALSE
    )
    fwrite(
        kmer_ref, 
        file = paste0("../../../GenomeAssembler_dev/data/", 
                      "QueryTable/QueryTable_kmer-", k,".csv"),
        showProgress = FALSE
    )    
}