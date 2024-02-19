suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults
pbapply::pboptions(char = "=")

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

# get org file and only import files that have single nucleotide breakpoints
org.file <- fread("../../data/org_file.csv")
org.file[, exp := paste0(`Fragmentation type`, "/", `Experiment folder`)]
org.file <- org.file[`RMSD?` == TRUE]

# org.file <- org.file[stringr::str_detect(string = org.file$exp, pattern = "TK6_Top2_mediated_DSBs/ETO_")]

# loop over outer breakage files
all.overlaps <- lapply(1:length(org.file$exp), function(outer.i){
    cat("Processing outer exp ", outer.i, "/", length(org.file$exp), ".\n", sep="")

    inner.ind <- 1:length(org.file$exp)
    inner.ind <- inner.ind[inner.ind != outer.i]

    # load in outer.i
    outer.file <- lapply(1:22, function(chr){
        load.outer.file <- paste0(
            "../../data/", 
            org.file$exp[outer.i], 
            "/breakpoint_positions/chr", chr, ".csv"
        )
        outer.file <- fread(load.outer.file, showProgress = FALSE)
        outer.file[, `:=`(seqnames = paste0("chr", chr), width = 1)]
        setcolorder(outer.file, c("seqnames", "start.pos", "width"))
        setnames(outer.file, c("seqnames", "start", "width"))
    })
    outer.file <- rbindlist(outer.file)
    outer.file <- plyranges::as_granges(outer.file)

    # nr of breakpoints
    outer.file.nr.breaks <- length(outer.file)

    # loop over inner breakage files
    inner.overlaps <- pbapply::pblapply(inner.ind, function(inner.i){
        # load in inner.i
        inner.file <- lapply(1:22, function(chr){
            load.inner.file <- paste0(
                "../../data/", 
                org.file$exp[inner.i], 
                "/breakpoint_positions/chr", chr, ".csv"
            )
            inner.file <- fread(load.inner.file, showProgress = FALSE)
            inner.file[, `:=`(seqnames = paste0("chr", chr), width = 1)]
            setcolorder(inner.file, c("seqnames", "start.pos", "width"))
            setnames(inner.file, c("seqnames", "start", "width"))
        })
        inner.file <- rbindlist(inner.file)
        inner.file <- plyranges::as_granges(inner.file)

        # nr of breakpoints
        inner.file.nr.breaks <- length(inner.file)

        # total number of breakpoints
        total_unique_breaks <- plyranges::bind_ranges(outer.file, inner.file)
        total_unique_breaks <- length(n_distinct(total_unique_breaks))

        # find overlaps
        overlap.files <- plyranges::filter_by_overlaps(outer.file, inner.file)
        overlaps <- length(overlap.files)/total_unique_breaks

        return(list(
            org.file$exp[inner.i], org.file$exp[outer.i], overlaps
        ))
    })

    gc(verbose = FALSE)
    inner.exp <- sapply(inner.overlaps, `[[`, 1)
    outer.exp <- sapply(inner.overlaps, `[[`, 2)
    overlaps <- sapply(inner.overlaps, `[[`, 3)

    return(list(
        inner.exp, outer.exp, overlaps
    ))
})

# combine results
df <- rbindlist(all.overlaps)
setnames(df, c("inner.exp", "outer.exp", "overlaps"))
dir.create(
    path = "../data/exact_breaks",
    showWarnings = FALSE
)
fwrite(df, "../data/exact_breaks/overlap_percent.csv")

# plot density distribution
p1 <- as_tibble(df) %>%
    dplyr::mutate(overlaps = overlaps*100) %>%
    ggplot(aes(x = overlaps)) + 
    geom_density(
        fill = "#69b3a2", 
        color = "black",
        alpha = 0.8) + 
    scale_x_continuous(limits = c(0, 100)) +
    theme_bw() + 
    theme(text = element_text(size = 20)) +
    coord_cartesian(xlim = c(0, 100)) +
    labs(
        x = "Overlaps, %",
        y = "Density",
        title = paste0("Percentage overlaps between identical \n", 
                       "breakpoint sites between all experiments")
    )

dir.create(
    path = "../figures/exact_breaks/",
    showWarnings = FALSE
)
ggsave(
    filename = "../figures/exact_breaks/BreakOverlaps.pdf",
    plot = p1,
    height = 6,
    width = 8
)