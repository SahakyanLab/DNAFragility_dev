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

        initialize = function(chr, which_exp_ind, control, seed){
            if(!missing(chr)) self$chr <- chr
            if(!missing(which_exp_ind)) self$which_exp_ind <- which_exp_ind
            if(!missing(control)) private$control <- control
            if(!missing(seed)) self$seed <- seed

            # get full org_file.csv
            private$get_org_file()
        },

        #' @description
        #' Calculate sequence-context effect via two adjacent
        #' base-pair of all aligned reads RMSD computations.
        #' @return None.
        calc_seq_effect = function(k = NULL, rmsd.range = c(-301, 301), 
                                   from_file = FALSE){
            rmsd.range <- rmsd.range[1]:rmsd.range[2]

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

                # load human reference genome
                private$get_ref(ind = i)

                # load breakpoint
                private$load_breakpoints(from_file = from_file)

                for(kmer in private$k){
                    cur.msg <- paste0("Calculating RMSD values for kmer ", kmer)
                    l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                    cat(cur.msg, l, "\n", sep = "")

                    # run RMSD calculations
                    private$run_rmsd(
                        rmsd.range = rmsd.range, 
                        k = kmer
                    )
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

        #' @field ref_rmsd character vector of a human reference genome chromosome.
        ref_rmsd = NULL,

        #' @field df_bp Data.Table of breakpoint positions for one chromosome.
        df_bp = NULL,

        #' @field k Numeric vector of k-mer size.
        k = NULL, 

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
        #' @return None.
        get_ref = function(ind){
            ref.seq <- private$org_file[ind, `Reference genome folder`]
            private$ref_rmsd <- paste0(
                "../../data/ref/", ref.seq, 
                "/chr", self$chr, ".fasta.gz"
            )
            private$ref <- Biostrings::readDNAStringSet(
                filepath = private$ref_rmsd
            )
        },

        #' @description 
        #' Load breakpoints for one chromosome.
        #' @return None.
        load_breakpoints = function(from_file = FALSE){
            fetch.file <- paste0(
                "../../data/", private$bp_exp,
                "/breakpoint_positions/chr", 
                self$chr, ".csv"
            )
            
            df <- fread(
                file = fetch.file,
                showProgress = FALSE
            )
            if("freq" %in% colnames(df)) df[, lev.dist := NULL]
            setorder(df, start.pos)

            if(private$control){
                # create control breakpoints
                if(from_file){
                    df <- fread(paste0(
                        "../lib/control_bp_kmer_8_seed_", self$seed, ".csv"
                    ))
                    df[, freq := 1]
                } else {
                    set.seed(self$seed)
                    to.sample <- 1500000
                    sample.points <- sample(
                        width(private$ref), 
                        size = to.sample, 
                        replace = FALSE
                    )
                    sample.points <- 
                        sample.points[!sample.points %in% df$start.pos]

                    # get substring of reference sequence
                    ref.seqs <- substring(
                        text = private$ref,
                        first = sample.points,
                        last = sample.points+1
                    )
                    to.keep <- which(!stringr::str_detect(
                        string = ref.seqs,
                        pattern = "N"
                    ))
                    sample.points <- sample.points[to.keep]
                    sample.points <- unique(sample.points)
                    sample.points <- sort(sample.points)
                    df <- data.table(
                        start.pos = sample.points, 
                        freq = 1
                    )
                }
            }
            private$df_bp <- df
        },

        #' @description
        #' Performs full process to obtain RMSD values.
        #' @param rmsd.range Numeric vector of range for RMSD calculations.
        #' @param kmer Numeric vector of kmer to perform calculation on.
        #' @return None.
        run_rmsd = function(rmsd.range, kmer){
            # calculate rmsd values 
            rmsd.values <- calc_kmer_freq(
                bp_pos = private$df_bp$start.pos,
                filename = private$ref_rmsd,
                k = kmer,
                rmsd_range = rmsd.range
            )

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
    )
)