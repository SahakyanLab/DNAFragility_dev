suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(rtracklayer)))
pbapply::pboptions(char = "=", type = "txt")
options(timeout=100000)

########################################################################
# load in all the motifs and positions from JASPAR
path_to_dat <- "../data/TFBS"

load_motif_mapping <- function(path){
    mapping <- readr::read_file(path) %>%
        str_split(., ">") %>%
        sapply(., function(x) {
            str_extract_all(x, "[.:A-Za-z0-9]*\\t[.:A-Za-z0-9]*(?=\\n)")
        }) %>%
        unlist() %>%
        sapply(., function(x) {
            str_split(x, "\\t")
        })
    names(mapping) <- sapply(mapping, function(x) {
        x[1]
    })
    mapping <- sapply(mapping, function(x) {
        x[2]
    })
    mapping
}

# Load mapping between matrix ID and motif name
file_to_download <- "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_redundant_pfms_jaspar.txt"
if(!file.exists("JASPAR2024_CORE_vertebrates_redundant_pfms_jaspar")){
    download.file(
        url = file_to_download,
        destfile = "JASPAR2024_CORE_vertebrates_redundant_pfms_jaspar.txt"
    )
}

jaspar_motifs <- load_motif_mapping("JASPAR2024_CORE_vertebrates_redundant_pfms_jaspar.txt")
hg38_chain <- import.chain("../../data/liftover/hg38ToHg19.over.chain")

# Load motifs bams into single object, and make method for anottations
if(!file.exists("jaspar_beds.bed")){
    download.file(
        url = "https://jaspar.elixir.no/download/data/2024/bed.tar.gz",
        destfile = "jaspar_beds"
    )
}
untar("jaspar_beds")
bed_path <- "./bed"

bed_files <- list.files(bed_path, ".bed", full.names = TRUE)
all_bed_names <- paste0(bed_path, "/", names(jaspar_motifs), ".bed")
valid_files <- intersect(all_bed_names, bed_files)
bed_ranges <- pbapply::pbsapply(valid_files, plyranges::read_bed)

# only keep hg19 and hg38 versions
assemblies <- lapply(1:length(bed_ranges), function(x){
    jaspar_label <- stringr::str_extract(
        names(bed_ranges[x]), 
        "bed/([A-Za-z0-9:.]*).bed", group = 1
    )
    tfbs_label <- as.character(jaspar_motifs[jaspar_label])

    # if any punctuations
    tfbs_label <- gsub("[[:punct:]]", "_", tfbs_label)
    tfbs_label <- gsub("___+", "", tfbs_label)
    tfbs_label <- gsub("__+", "_", tfbs_label)

    # capitalise
    tfbs_label <- toupper(tfbs_label)

    data.table(
        ind = x,
        assembly = stringr::str_extract(mcols(bed_ranges[[x]])$name[[1]], ".+?(?=_)"),
        jaspar_name = jaspar_label,
        tfbs_label = tfbs_label
    )
})
assemblies <- as.data.table(do.call(rbind, assemblies))

# if multiple TFBS exist, keep the first occurring one.
assemblies <- as_tibble(assemblies) %>% 
    dplyr::filter(grepl("hg19|hg38", assembly, ignore.case = TRUE)) %>% 
    dplyr::arrange(tfbs_label, jaspar_name) %>% 
    dplyr::group_by(tfbs_label) %>% 
    dplyr::slice(1) %>% 
    dplyr::ungroup() %>% 
    as.data.frame()

dummy.org.file <- data.table(
    `Fragmentation type` = "02_JASPAR",
    `Experiment folder` = rep("", nrow(assemblies)),
    `Reference genome folder` = "",
    Processed = TRUE,
    Index = 0,
    Start = 0,
    End = 0,
    `Both Ends?` = TRUE,
    `File format` = "CSV",
    Category_Main = "",
    Category_general = "TFBS_JASPAR",
    Category_Hex = "#619B61",
    Class = "TFBS",
    `DSB Map` = FALSE,
    Alignment_strand = "both",
    `RMSD?` = FALSE,
    kmer_2 = "", 
    kmer_4 = "",
    kmer_6 = "",
    kmer_8 = ""
)

# liftover to hg19
for(x in 1:nrow(assemblies)){
    cat(paste0("Processing TF ", x, "/", nrow(assemblies), ".\n"))

    genome_assembly <- stringr::str_extract(mcols(bed_ranges[[assemblies$ind[x]]])$name[[1]], ".+?(?=_)")
    bed_grange <- bed_ranges[[assemblies$ind[x]]]
    bed_grange_label <- assemblies$tfbs_label[x]

    if(grepl("hg38", genome_assembly)){
        liftover_hg38 <- suppressMessages(liftOver(bed_grange, hg38_chain))
        bed_grange <- unlist(as(liftover_hg38, "GRangesList"))
    } 

    liftover_hg38 <- bed_grange %>% 
        dplyr::select(-name, -score) %>% 
        plyranges::filter(seqnames %in% paste0("chr", 1:22)) %>%
        reduce(.)

    dt_this_tfbs <- as.data.table(liftover_hg38)
    dt_this_tfbs[, `:=`(
        seqnames = as.character(seqnames), 
        width = NULL,
        strand = ifelse(strand == "-", "-", "+")
    )]
    setnames(dt_this_tfbs, c("Chromosome", "Start", "End", "Strand"))

    # save file for use in kmertone
    dir.create(
        path = paste0("./", bed_grange_label, "/kmertone"),
        showWarnings = FALSE,
        recursive = TRUE
    )

    # save chromosome separated files
    for(chr in unique(dt_this_tfbs$Chromosome)){
        fwrite(
        dt_this_tfbs[Chromosome == chr],
        file = paste0("./", bed_grange_label, "/kmertone/", chr, ".csv")
        )
    }

    fwrite(
        dt_this_tfbs,
        file = paste0("./", bed_grange_label, "/", bed_grange_label, ".csv")
    )

    # update dummy org file
    dummy.org.file$`Experiment folder`[x] <- bed_grange_label
    dummy.org.file$`Reference genome folder`[x] <- genome_assembly
    dummy.org.file$Category_Main[x] <- paste0("TFBS_", bed_grange_label)
}

fwrite(
    dummy.org.file,
    "../dummy_org_file.csv"
)

# # rbind dummy matrix to org.file but only do the below if you're 100% sure it works
# org.file <- rbind(org.file, dummy.org.file)
# fwrite(
#     org.file,
#     "../org_file.csv"
# )