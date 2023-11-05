SequenceEffect <- R6::R6Class(
    classname = "SequenceEffect",
    public = list(
        #' @field chr Numeric vector of chromosome number.
        chr = 1,

        #' @field which_exp_ind Numeric vector of the index of the experiment 
        #' to process from org_file. If NULL, all exp will be processed.
        which_exp_ind = 1,

        #' @field plots list of ggplot objects showing the RMSD values vs. 
        #' positions away from breakpoint per kmer.
        plots = NULL,

        #' @field seed Numeric vector to fix the seed.
        seed = 1234,

        #' @field results List of results stored if return_vals is TRUE.
        results = NULL,

        initialize = function(chr, which_exp_ind, control, seed, return_vals){
            if(!missing(chr)) self$chr <- chr
            if(!missing(which_exp_ind)) self$which_exp_ind <- which_exp_ind
            if(!missing(control)) private$control <- control
            if(!missing(seed)) self$seed <- seed
            if(!missing(return_vals)) private$return_vals <- return_vals

            if(private$return_vals){
                self$results <- vector(
                    mode = "list", 
                    length = length(self$chr)
                )
                names(self$results) <- paste0("chr", self$chr)
            }

            # get full org_file.csv
            private$get_org_file()
        },

        #' @description
        #' Calculate sequence-context effect via two adjacent
        #' base-pair of all aligned reads RMSD computations.
        #' @return None.
        calc_seq_effect = function(k = NULL, rmsd.range = c(-301, 301), 
                                   from_file = FALSE, break_type = NULL){
            private$rmsd_range <- rmsd.range[1]:rmsd.range[2]
            if(!is.null(break_type)) private$break_type <- break_type

            if(is.null(self$which_exp_ind)){
                self$which_exp_ind <- which(
                    (private$org_file$`DSB Map` == "TRUE") & 
                    (private$org_file$`RMSD?` == "TRUE")
                )
            }
            len.of.loop <- self$which_exp_ind

            # if k is not specified, perform calculation over c(4,6,8).
            if(length(k) == 1){
                private$k <- k
            } else {
                private$k <- c(2, 4, 6, 8)
            }

            # loop over each experiment
            for(i in len.of.loop){
                start.time <- Sys.time()

                private$bp_exp <- paste0(
                    private$org_file[i, `Fragmentation type`], "/",
                    private$org_file[i, `Experiment folder`]
                )

                cur.msg <- paste0("Processing exp ", i, "/", 
                                  nrow(private$org_file),
                                  ". ", private$bp_exp)
                cat(cur.msg, "\n", sep = "")

                for(chr in self$chr){
                    t1 <- Sys.time()
                    cur.msg <- paste0("Calculating RMSD values for chromosome: ", chr)
                    l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                    cat(paste0(cur.msg, l))

                    # load human reference genome
                    private$get_ref(ind = i, chr = chr)

                    # load breakpoint
                    private$load_breakpoints(chr = chr, from_file = from_file)

                    for(kmer in private$k){
                        if(length(private$k) > 1){
                            cur.msg <- paste0("Calculating RMSD values for kmer ", kmer)
                            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                            cat(cur.msg, l, "\n", sep = "")
                        }

                        # run RMSD calculations
                        private$run_rmsd(k = kmer, chr = chr)
                    }

                    if(length(private$k) <= 1){
                        total.time <- Sys.time() - t1
                        cat("DONE! --", signif(total.time[[1]], 2), 
                            attr(total.time, "units"), "\n")
                    }
                }

                # time taken for full processing for this experiment
                final.t <- Sys.time() - start.time
                cat(paste(c(rep("-", 70), "\n"), collapse = ""))
                cat("Final time taken:", signif(final.t[[1]], digits = 2), 
                    attr(final.t, "units"), "\n")
                cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            }
        },

        #' @description
        #' Generates ggplots of calculated RMSD values.
        #' @return None.
        quick_plot_check = function(){
            t1 <- Sys.time()
            cur.msg <- "Generating RMSD plots"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            plots.saved <- c()
            for(kmer in private$k){
                files <- list.files(
                    path = paste0("../data/", private$bp_exp),
                    pattern = paste0(
                        ifelse(private$control,
                        paste0("^control_rmsd_kmer_", kmer, 
                               "_seed_", self$seed), 
                        paste0("^rmsd_kmer_", kmer))
                    ),
                    full.names = TRUE
                )
                df <- readRDS(file = files)
                limits <- length(df)/2-1
                df <- as_tibble(df) %>% 
                    dplyr::mutate( 
                        x = -limits:(length(df)-limits-1)) %>% 
                    dplyr::rename(y = value)

                p1 <- df %>%
                    ggplot(aes(x = x, y = y)) + 
                    geom_line(linewidth = 0.8) + 
                    theme_bw() + 
                    labs(
                        x = "Position away from breakpoint",
                        y = "RMSD",
                        title = paste0("Exp: ", private$bp_exp),
                        subtitle = paste0("kmer_", kmer)
                    ) +
                    ylim(0, NA)

                self$plots[[paste0("kmer_", kmer)]] <- p1
                plots.saved <- c(plots.saved, paste0("self$plots$kmer_", kmer))
            }

            # time taken for full processing for this experiment
            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], digits = 2), 
                attr(total.time, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cat("View plots:", plots.saved, "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))            
        }
    ), 
    private = list(
        #' @field control Boolean. If TRUE, will perform RMSD on controls.
        control = FALSE,

        #' @field org_file Data.table of the full org_file.csv
        org_file = NULL,

        #' @field bp_exp Character vector of fragment type and experiment.
        bp_exp = NULL,

        #' @field ref DNAStringSet of a human reference genome chromosome.
        ref = NULL,

        #' @field chrs_len Numeric vector of the length of each chromosome of the reference sequence.
        chrs_len = NULL,

        #' @field ref_rmsd character vector of a human reference genome chromosome.
        ref_rmsd = NULL,

        #' @field df_bp Data.Table of breakpoint positions for one chromosome.
        df_bp = NULL,

        #' @field k Numeric vector of k-mer size.
        k = NULL, 

        #' @field rmsd.range Numeric vector of range for RMSD calculations.
        rmsd_range = NULL,

        #' @field return_vals Boolean. If TRUE, will return values instead of saving as file.
        return_vals = FALSE,

        #' @field break_type Character vector. Choose the long range cut-off for the breakage class. 
        break_type = NULL,

        #' @description
        #' Import full org_file.csv and filter for rows to be processed.
        #' @return None.
        get_org_file = function(){
            private$org_file <- fread(
                "../../data/org_file.csv",
                showProgress = FALSE
            )            
        },

        #' @description
        #' Get human reference genome file.
        #' @param ind Numeric vector of index subsetting org_file.
        #' @param chr Numeric vector of chromosome number.
        #' @return None.
        get_ref = function(ind, chr){
            ref.seq <- private$org_file[ind, `Reference genome folder`]
            private$ref_rmsd <- paste0(
                "../../data/ref/", ref.seq, 
                "/chr", chr, ".fasta.gz"
            )
            
            if(private$control){
                private$ref <- Biostrings::readDNAStringSet(
                    filepath = private$ref_rmsd
                )
            }

            # get the start and end positions of each chromosome
            temp <- switch(ref.seq,
                "hg18" = BSgenome.Hsapiens.UCSC.hg18::BSgenome.Hsapiens.UCSC.hg18,
                "hg19" = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                "hg38" = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                "hs37d5" = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
            )
            refseq.table <- as.data.frame(temp@seqinfo)
            refseq.table <- refseq.table[grepl(
                pattern = "^chr([1-9]|1[0-9]|2[0-2])$", 
                x = rownames(refseq.table)
            ),]

            if(nrow(refseq.table) == 0){
                refseq.table <- as.data.frame(temp@seqinfo)
                refseq.table <- refseq.table[grepl(
                    pattern = "^([1-9]|1[0-9]|2[0-2])$", 
                    x = rownames(refseq.table)
                ),]
                rownames(refseq.table) <- paste0("chr", rownames(refseq.table))
            }

            private$chrs_len <- refseq.table$seqlengths         
            private$chrs_len <- setNames(private$chrs_len, rownames(refseq.table))
        },

        #' @description 
        #' Load breakpoints for one chromosome.
        #' @param from_file Boolean. If True, will use control breakpoints obtained via Kmertone.
        #' @param chr Numeric vector of chromosome number.
        #' @return None.
        load_breakpoints = function(from_file = FALSE, chr){
            fetch.file <- paste0(
                "../../data/", private$bp_exp,
                "/breakpoint_positions/chr", 
                chr, ".csv"
            )
            
            df <- fread(
                file = fetch.file,
                showProgress = FALSE
            )
            if("freq" %in% colnames(df)) df[, lev.dist := NULL]
            setorder(df, start.pos)

            # make sure no position will be out-of-bounds
            df[, `:=`(
                lower_end = start.pos+min(private$rmsd_range),
                upper_end = start.pos+max(private$rmsd_range)
            )]
            df[, to_keep := ifelse(
                lower_end > 0 & upper_end <= private$chrs_len[chr],
                TRUE, FALSE
            )]
            df <- df[to_keep == TRUE]
            df[, `:=`(lower_end = NULL, upper_end = NULL, to_keep = NULL)]

            if(private$control){
                round_to_nearest_even <- function(x) round(x/2)*2
                ranges <- fread(paste0(
                    "../../../04_DNAFragility/data/range_effects/", 
                    "MaxValuesFromClustersByType.csv"
                ))
                row_id <- which(ranges$type == private$break_type)
                all_ranges <- ranges[row_id,-"type"]
                all_ranges <- as.list(all_ranges)

                # round to the nearest even number
                long_range <- round_to_nearest_even(all_ranges$mid.range)

                if(from_file){
                    # how many positions to sample from (max)
                    to_sample_max <- nrow(df)*0.4
                    if(grepl("Ultrasonication", private$bp_exp) | 
                        (private$break_type == "biological")
                        ){
                        to_sample_max <- 5000000
                    }

                    # from control regions
                    fetch.file <- paste0(
                        "../../04_KmericAnalysis/data/kmertone/", 
                        private$bp_exp,
                        "/control_coordinates/chr", 
                        chr, ".csv"
                    )

                    sample_points <- fread(
                        file = fetch.file,
                        showProgress = FALSE
                    )
                    sample_points[, `:=`(
                        seqnames = paste0("chr", chr),
                        width = end-start
                    )]
                    setcolorder(sample_points, c("seqnames", "start", "end", "width"))

                    # sample rows proportional to width
                    set.seed(self$seed)
                    sample_rows <- sample(
                        nrow(sample_points), 
                        size = to_sample_max, 
                        replace = TRUE, 
                        prob = sample_points$width
                    )
                    sample_points <- sample_points[sample_rows]

                    # sample within each chosen range
                    sample_points[, position := sample(start:end, size = 1, replace = FALSE), 
                                by = 1:dim(sample_points)[1]]
                    sample_points[, start := position]
                    sample_points[, c("seqnames", "end", "width", "position") := NULL]
                    setorder(sample_points, start)
                    setnames(sample_points, "start.pos")

                    # remove points that may overlap within 500 bases from true break
                    sample_points[, `:=`(
                        seqnames = paste0("chr", chr),
                        start = start.pos-long_range,
                        end = start.pos+long_range
                        # start = start.pos-ceiling(length(private$rmsd_range)/2-1),
                        # end = start.pos+ceiling(length(private$rmsd_range)/2-1)
                    )]
                    sample_points_granges <- plyranges::as_granges(sample_points)
                    
                    setnames(df, "start")
                    df[, `:=`(
                        seqnames = paste0("chr", chr),
                        width = 1
                    )]
                    df <- plyranges::as_granges(df)

                    sample_points <- plyranges::filter_by_non_overlaps(
                        sample_points_granges, 
                        df
                    )
                    df <- as.data.table(mcols(sample_points)$start.pos)
                    setnames(df, "start.pos")
                } else {
                    # how many positions to sample from (max)
                    to_sample_max <- 20000000

                    # random sampling
                    set.seed(self$seed)
                    sample_points <- sample(
                        width(private$ref),
                        size = to_sample_max, 
                        replace = FALSE
                    )

                    # remove points that may overlap within 500 bases from the true break
                    sample_points_granges <- data.table(
                        seqnames = paste0("chr", chr),
                        start = sample_points-long_range,
                        end = sample_points+long_range,
                        # start = sample_points-ceiling(length(private$rmsd_range)/2-1),
                        # end = sample_points+ceiling(length(private$rmsd_range)/2-1), 
                        start.pos = sample_points
                    ) %>% plyranges::as_granges()
                    
                    setnames(df, "start")
                    df[, `:=`(
                        seqnames = paste0("chr", chr),
                        width = 1
                    )]
                    df <- plyranges::as_granges(df)

                    sample_points <- plyranges::filter_by_non_overlaps(
                        sample_points_granges, 
                        df
                    )
                    sample_points <- as.numeric(mcols(sample_points)$start.pos)

                    # get substring of reference sequence
                    ref_seqs <- substring(
                        text = private$ref,
                        first = sample_points,
                        last = sample_points+1
                    )
                    to.keep <- which(!stringr::str_detect(
                        string = ref_seqs,
                        pattern = "N"
                    ))
                    sample_points <- sample_points[to.keep]
                    sample_points <- unique(sample_points)
                    sample_points <- sort(sample_points)
                    df <- data.table(start.pos = sample_points)
                }
                df <- unique(df, by = "start.pos")
            } 
            private$df_bp <- df
        },

        #' @description
        #' Performs full process to obtain RMSD values.
        #' @param kmer Numeric vector of kmer to perform calculation on.
        #' @param chr Numeric vector of chromosome number.
        #' @return None.
        run_rmsd = function(kmer, chr){
            # calculate rmsd values 
            rmsd.values <- calc_kmer_freq(
                bp_pos = private$df_bp$start.pos,
                filename = private$ref_rmsd,
                kmer = kmer,
                rmsd_range = private$rmsd_range
            )

            if(private$return_vals){
                self$results[[paste0("chr", chr)]] <- rmsd.values
            } else {
                # save absolute frequency values for motif analysis
                dir.create(
                    path = paste0("../data/", private$bp_exp),
                    showWarnings = FALSE,
                    recursive = TRUE
                )
                saveRDS(
                    object = rmsd.values,
                    file = paste0(
                        "../data/", private$bp_exp,
                        ifelse(private$control, 
                        paste0("/control_rmsd_kmer_", 
                            kmer, "_seed_", self$seed), 
                        paste0("/rmsd_kmer_", kmer)),
                        ".Rdata")
                )
            }
        }
    )
)