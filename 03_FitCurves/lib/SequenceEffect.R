SequenceEffect <- R6::R6Class(
    classname = "SequenceEffect",
    public = list(
        #' @field chr Numeric vector of chromosome number.
        chr = NULL,

        #' @field which_exp_ind Numeric vector of the index of the experiment 
        #' to process from org_file. If NULL, all exp will be processed.
        which_exp_ind = NULL,

        #' @field cores Numeric vector of the nr of CPUs to use.
        cores = NULL,

        initialize = function(chr, which_exp_ind, cores, control){
            if(!missing(chr)) self$chr <- chr
            if(!missing(which_exp_ind)) self$which_exp_ind <- which_exp_ind
            if(!missing(cores)) self$cores <- cores
            if(!missing(control)) private$control <- control

            # get full org_file.csv
            private$get_org_file()
        },

        #' @description
        #' Calculate sequence-context effect via two adjacent
        #' base-pair of all aligned reads RMSD computations.
        #' @return None.
        calc_seq_effect = function(k, rmsd.range = c(-301, 301)){
            rmsd.range <- rmsd.range[1]:rmsd.range[2]

            if(is.null(self$which_exp_ind)){
                len.of.loop <- 1:nrow(private$org_file)
            } else {
                len.of.loop <- self$which_exp_ind
            }

            # if k is not specified, perform calculation over c(2,4,6,8).
            if(length(k) == 1){
                private$k <- k
            } else {
                private$k <- c(2, 4, 6, 8)
            }

            # loop over each experiment
            for(i in len.of.loop){
                start.time <- Sys.time()
                cur.msg <- paste0("Processing experiment ", i, 
                                  " of ", nrow(private$org_file))
                l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                cat(cur.msg, l, "\n", sep = "")

                private$bp_exp <- paste0(
                    private$org_file[i, `Fragmentation type`], "/",
                    private$org_file[i, `Experiment folder`]
                )

                # load breakpoint
                private$load_breakpoints()

                # load human reference genome
                private$get_ref(ind = i)

                for(kmer in private$k){
                    cur.msg <- paste0("Calculating RMSD values for kmer ", kmer)
                    l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                    cat(cur.msg, l, "\n", sep = "")

                    # generate k-mer table
                    private$generate_table(kmer = kmer)

                    # run RMSD calculations
                    private$run_rmsd(
                        rmsd.range = rmsd.range, 
                        k = kmer
                    )
                }

                # time taken for full processing for this experiment
                final.t <- Sys.time() - start.time
                cat(paste(c(rep("-", 70), "\n"), collapse = ""))
                cat("Final time taken:", signif(final.t[[1]], digits = 3), 
                    attr(final.t, "units"), "\n")
                cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            }
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

        #' @field df_bp Data.Table of breakpoint positions for one chromosome.
        df_bp = NULL,

        #' @field kmer_list Character vector of k-mers.
        kmer_list = NULL,

        #' @field kmer_ref Data.table of first occurring 
        #' kmer in lexicological order.
        kmer_ref = NULL,

        #' @field k Numeric vector of k-mer size.
        k = NULL,

        #' @description
        #' Import full org_file.csv and filter for rows to be processed.
        #' @return None.
        get_org_file = function(){
            df <- fread(
                "../../data/org_file.csv",
                showProgress = FALSE
            )
            private$org_file <- df[`DSB Map` == "TRUE"]
        },

        #' @description 
        #' A utility function to generate k-mers.
        #' Then, only keeps first occurring k-mer in lexicological order.
        #' @return None.
        generate_table = function(kmer){
            k.mers <- do.call(data.table::CJ, 
                              rep(list(c("A", "C", "G", "T")), kmer))
            private$kmer_list <- k.mers[, do.call(paste0, .SD)]
            rev.comp <- as.character(
                Biostrings::reverseComplement(Biostrings::DNAStringSet(private$kmer_list))
            )
            kmer_ref <- data.table('fwd' = private$kmer_list, 'rev.comp' = rev.comp)
            kmer_ref[, cond := ifelse(seq(1:nrow(.SD)) < match(fwd, rev.comp), 
            TRUE, ifelse(fwd == rev.comp, TRUE, FALSE))]
            private$kmer_ref <- kmer_ref[cond == TRUE, .(fwd, rev.comp)]
        },

        #' @description
        #' Get human reference genome file.
        #' @param ind Numeric vector of index subsetting org_file.
        #' @return None.
        get_ref = function(ind){
            ref.seq <- private$org_file[ind, `Reference genome folder`]
            private$ref <- Biostrings::readDNAStringSet(
                filepath = paste0("../../data/ref/", ref.seq, 
                                  "/chr", self$chr, ".fasta.gz")
            )
        },

        #' @description 
        #' Load breakpoints for one chromosome.
        #' @return None.
        load_breakpoints = function(){
            fetch.file <- paste0(
                "../../data/", private$bp_exp,
                "/breakpoint_positions/chr", 
                self$chr, ".csv"
            )
            
            df <- fread(
                file = fetch.file,
                select = c("start.pos", "freq"),
                showProgress = FALSE
            )
            setorder(df, start.pos)

            if(private$control){
                # create control breakpoints
                set.seed(seed = 1234)
                df <- data.table(
                    start.pos = floor(runif(
                        n = nrow(df), 
                        min = 1, 
                        max = width(private$ref))))[order(start.pos)]
            }
            private$df_bp <- df
        },

        #' @description
        #' Performs full process to obtain RMSD values.
        #' @param rmsd.range Numeric vector of range for RMSD calculations.
        #' @param kmer Numeric vector of kmer to perform calculation on.
        #' @return None.
        run_rmsd = function(rmsd.range, kmer){
            # kmer frequency calculations
            freq <- pbapply::pblapply(rmsd.range, function(x){
                freq.vals <- private$calc_kmer_freq(ind = x, kmer = kmer)
                norm.vals <- freq.vals/sum(freq.vals, na.rm = TRUE)
                return(list(freq.vals, norm.vals))
            }, cl = self$cores)
            freq.vals <- sapply(freq, `[[`, 1)
            norm.vals <- sapply(freq, `[[`, 2)

            # save absolute frequency values for motif analysis
            dir.create(
                path = paste0("../data/", private$bp_exp),
                showWarnings = FALSE,
                recursive = TRUE
            )
            saveRDS(
                object = freq.vals,
                file = paste0("../data/", private$bp_exp,
                              ifelse(private$control, 
                              "/control_freq_rmsd_kmer_", 
                              "/freq_rmsd_kmer_"),
                              kmer, ".Rdata")
            )
            # save normalised frequency values
            saveRDS(
                object = norm.vals,
                file = paste0("../data/", private$bp_exp,
                              ifelse(private$control, 
                              "/control_rmsd_kmer_", 
                              "/normalised_freq_rmsd_kmer_"),
                              kmer, ".Rdata")
            )

            # save rmsd values
            rmsd.values <- pbapply::pbsapply(1:(dim(norm.vals)[2]-1), function(x){
                private$calc_rmsd(norm.vals[, x], norm.vals[, (x+1)])
            })
            saveRDS(
                object = rmsd.values,
                file = paste0("../data/", private$bp_exp,
                              ifelse(private$control, 
                              "/control_rmsd_kmer_", 
                              "/rmsd_kmer_"),
                              kmer, ".Rdata")
            )
        },

        #' @description
        #' Calculate kmer frequencies at each shifted bp position.
        #' @param ind Numeric vector of index subsetting org_file.
        #' @param kmer Numeric vector of kmer to perform calculation on.
        #' @return Data.Table of kmer frequencies.
        calc_kmer_freq = function(ind, kmer){
            dfcopy <- copy(private$df_bp)
            dfcopy[, `:=`(start.pos = start.pos + ind)]
            
            if((kmer %% 2) == 0){
                ending.pos <- ceiling((kmer-1)/2)
                starting.pos <- (kmer-1)-ending.pos
            } else {
                interval <- (kmer-1)/2
                starting.pos <- ending.pos <- interval
            }

            # obtain k-mers from alignment data
            dfcopy[,`:=`(start = start.pos - starting.pos, end = start.pos + ending.pos)]
            dfcopy[, start.pos := NULL]

            # extract k-meric counts and relative frequencies
            dfcopy[, `:=`(fwd = stringr::str_sub(
                string = private$ref, start = start, end = end))]

            if("-" %in% unique(dfcopy$strand)){
                dfcopy[strand == "-", fwd := paste(
                    Biostrings::reverseComplement(Biostrings::DNAStringSet(fwd))
                )]
            }
            dfcopy <- dfcopy[, .(n = .N), by = .(fwd)][!stringr::str_detect(
                string = fwd, pattern = "N")]
            
            # account for strand symmetry
            fwd.ind <- match(private$kmer_ref$fwd, dfcopy$fwd)
            rev.comp.ind <- match(private$kmer_ref$rev.comp, dfcopy$fwd)
            
            # update data frame with dyad frequency count
            private$kmer_ref[, `:=`(freq = dfcopy$n[fwd.ind] + dfcopy$n[rev.comp.ind])]
            res <- private$kmer_ref$freq
            private$kmer_ref[, freq := NULL]

            # return data frame of k-mers with associated breakpoint frequencies
            return(res)
        },

        #' @description
        #' Calculate root mean squared deviation between two points.
        #' @param a Numeric vector of k-mer frequencies.
        #' @param b Numeric vector of k-mer frequencies.
        #' @return Numeric vector of RMSD value.
        calc_rmsd = function(a, b){
            return(sqrt(mean((a-b)^2, na.rm = TRUE)))
        }
    )
)