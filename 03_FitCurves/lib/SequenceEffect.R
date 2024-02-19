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

        size_of_data = NULL,

        initialize = function(chr, which_exp_ind, control, seed, 
                              return_vals, assembly){
            if(!missing(chr)) self$chr <- chr
            if(!missing(which_exp_ind)) self$which_exp_ind <- which_exp_ind
            if(!missing(control)) private$control <- control
            if(!missing(seed)) self$seed <- seed
            if(!missing(return_vals)) private$return_vals <- return_vals

            if(private$return_vals){
                self$size_of_data <- self$results <- vector(
                    mode = "list", 
                    length = length(self$chr)
                )
                names(self$size_of_data) <- names(self$results) <- paste0("chr", self$chr)
            }

            if(is.null(self$which_exp_ind) & private$control){
                private$pure_control_study <- TRUE
            }

            # get full org_file.csv
            private$get_org_file()
        },

        #' @description
        #' Calculate sequence-context effect via two adjacent
        #' base-pair of all aligned reads RMSD computations.
        #' @return None.
        calc_seq_effect = function(k = NULL, rmsd.range = c(-301, 301), 
                                   from_file = FALSE, break_type = NULL,
                                   break_data = NULL, assembly = NULL){
            private$rmsd_range <- rmsd.range[1]:rmsd.range[2]
            if(!is.null(break_type)) private$break_type <- break_type

            if(!is.null(break_data)){
                for(i in 1:length(break_data)){
                    private$break_data[[names(break_data)[i]]] <- break_data[[names(break_data)[i]]]
                }
            }

            if(is.null(private$break_data)){
                if(is.null(self$which_exp_ind)){
                    self$which_exp_ind <- which(
                        (private$org_file$`DSB Map` == "TRUE") & 
                        (private$org_file$`RMSD?` == "TRUE")
                    )
                }
                len.of.loop <- self$which_exp_ind

                if(private$pure_control_study) len.of.loop <- 1
            } else {
                len.of.loop <- 1
            }

            # if k is not specified, perform calculation over c(4,6,8).
            if(length(k) == 1){
                private$k <- k
            } else {
                private$k <- c(2, 4, 6, 8)
            }

            # loop over each experiment
            for(i in len.of.loop){
                start.time <- Sys.time()

                if(is.null(private$break_data)){
                    if(private$pure_control_study){
                        private$assembly <- "hg19"
                        cur.msg <- paste0("Running a pure control study.")
                    } else {
                        private$bp_exp <- paste0(
                            private$org_file[i, `Fragmentation type`], "/",
                            private$org_file[i, `Experiment folder`]
                        )
                        private$assembly <- private$org_file[i, `Reference genome folder`]
                        
                        cur.msg <- paste0("Processing ", match(i, len.of.loop), "/", 
                                        length(len.of.loop),
                                        ". ", private$bp_exp)
                    }
                } else {
                    private$bp_exp <- NULL
                    private$assembly <- ifelse(is.null(assembly), "hg19", assembly)

                    cur.msg <- "Processing current experiment."
                }
                cat(cur.msg, "\n", sep = "")

                for(chr in self$chr){
                    t1 <- Sys.time()
                    cur.msg <- paste0("Calculating RMSD values for chromosome: ", chr)
                    l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                    cat(paste0(cur.msg, l))

                    # load human reference genome
                    private$get_ref(ind = i, chr = chr)

                    for(kmer in private$k){
                        if(length(private$k) > 1){
                            cur.msg <- paste0("Calculating RMSD values for kmer ", kmer)
                            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                            cat(cur.msg, l, "\n", sep = "")
                        }

                        # load breakpoint
                        private$load_breakpoints(
                            chr = chr, 
                            from_file = from_file,
                            k = kmer
                        )

                        # run RMSD calculations
                        private$run_rmsd(kmer = kmer, chr = chr)
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

        #' @field break_data Numeric vector of RMSD values. 
        break_data = NULL,

        #' @field assembly Character vector. Select version of the human genome assembly.
        assembly = "hg19",

        pure_control_study = FALSE,

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
            if(is.null(private$break_data)){
                ref.seq <- private$org_file[ind, `Reference genome folder`]
            } else {
                ref.seq <- private$assembly
            }
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
        #' @param kmer Numeric vector of kmer to perform calculation on.
        #' @return None.
        load_breakpoints = function(chr, kmer, from_file = FALSE){
            to_sample_max <- 30000000
            half_width <- max(abs(max(private$rmsd_range)))+ceiling(kmer/2)

            if(!private$pure_control_study){
                if(is.null(private$break_data)){
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
                } else {
                    df <- private$break_data[[paste0("chr", chr)]]
                    setnames(df, "start.pos")
                }
                setorder(df, start.pos)

                # make sure no position will be out-of-bounds
                df[, `:=`(
                    lower_end = start.pos-half_width,
                    upper_end = start.pos+half_width
                )]
                df[, to_keep := ifelse(
                    lower_end > 0 & upper_end < private$chrs_len[chr],
                    TRUE, FALSE
                )]
                df <- df[to_keep == TRUE]
                df[, `:=`(lower_end = NULL, upper_end = NULL, to_keep = NULL)]
                
                # how many positions to sample from (max) if 
                # performing a comparative control study
                to_sample_max <- nrow(df)
            }

            if(private$control){
                if(from_file){                    
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
                    setcolorder(sample_points, c("start", "end"))
                    sample_points[, seqnames := paste0("chr", chr)]
                    sample_points[, chr_len := private$chrs_len[seqnames]]
                    
                    sample_points[, `:=`(
                        start = pmax(start, half_width),
                        end = pmin(end, chr_len-half_width)
                    )]
                    sample_points[, width := end-start]
                    sample_points[, to_keep := ifelse(width > 0, TRUE, FALSE)]
                    sample_points <- sample_points[to_keep == TRUE]
                    sample_points[, to_keep := NULL]

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
                    if(length(chr) == 1){
                        sample_points[, seq_start_ID := start]
                    } else {
                        sample_points[, seq_start_ID := paste0(seqnames, "_", start)]
                    }

                    # only keep the first N elements if a group has more entries than there 
                    # are positions to be sampled from within a control region. 
                    # sample_points[, position := start+sample(width[1], size = 1, replace = FALSE), by = seq_start_ID]
                    # df <- sample_points
                    
                    df <- sample_points[, 
                        .SD[sample(.N, pmin(.N, width))], 
                        by = seq_start_ID
                    ]

                    df[, position := start+sample(
                        width[1], size = .N, replace = FALSE
                    ), by = seq_start_ID]
                    df[, start := position]
                    df[, c(
                        "seqnames", 
                        "end", 
                        "width", 
                        "position",
                        "chr_len", 
                        "seq_start_ID"
                    ) := NULL]
                    setorder(df, start)
                    setnames(df, "start.pos")
                    setorder(df, start.pos)
                } else {
                    # random sampling of regions not out-of-bounds
                    # not overlapping with true breaks
                    sample_points <- pot_range <- GenomicRanges::GRanges(
                        seqnames = paste0("chr", chr),
                        IRanges(
                            start = (half_width):(width(private$ref)-half_width),
                            width = 1
                        )
                    )
                    
                    if(!private$pure_control_study){
                        range_to_avoid <- GenomicRanges::GRanges(
                            seqnames = paste0("chr", chr),
                            IRanges(
                                start = df$start.pos-half_width,
                                end = df$start.pos+half_width
                            )
                        )

                        sample_points <- plyranges::filter_by_non_overlaps(
                            pot_range, range_to_avoid
                        )
                    }
                    pot_range <- start(sample_points)

                    set.seed(self$seed)
                    sample_points <- sample(
                        x = pot_range,
                        size = min(to_sample_max, length(pot_range)), 
                        replace = FALSE
                    )
                    df <- data.table(start.pos = sample_points)
                    setorder(df, start.pos)
                }
                df <- unique(df, by = "start.pos")
            } 
            private$df_bp <- df
            self$size_of_data[[paste0("chr", chr)]] <- nrow(private$df_bp)
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

            # if performing both a true and control study, then normalise
            # by the bigger denominator {true, control}_breakpoint data set
            # to compare both on the same set of RMSD y-axis
            # if(private$control){
            #     rmsd.values <- rmsd.values/self$size_of_data
            # }

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
                        "/chr", chr,
                        ifelse(private$control, 
                        paste0("_control_rmsd_kmer_", kmer, "_seed_", self$seed), 
                        paste0("_rmsd_kmer_", kmer)),
                        ".Rdata"
                    )
                )
            }
        }
    )
)