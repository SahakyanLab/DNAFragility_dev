AnalyseKmers <- R6::R6Class(
    classname = "AnalyseKmers",
    public = list(
        #' @field which_exp_ind Numeric vector of the index of the experiment 
        #' to process from org_file. If NULL, all exp will be processed.
        which_exp_ind = NULL,

        #' @field cores Numeric vector of the nr of CPUs to use.
        cores = NULL,

        #' @field kmer_list Character vector of k-mers.
        kmer_list = NULL,

        #' @field kmer_ref Data.table of first occurring 
        #' kmer in lexicological order.
        kmer_ref = NULL,

        initialize = function(which_exp_ind, cores){
            if(!missing(which_exp_ind)) self$which_exp_ind <- which_exp_ind
            if(!missing(cores)) self$cores <- cores

            # get full org_file.csv
            private$get_org_file()
        },

        #' @description
        #' Run kmertone software and obtain weight factors for each kmer.
        #' @param k Numeric vector of kmer to perform calculation on.
        #' @return None.
        get_kmer_enrichment = function(k = NULL){
            if(is.null(self$which_exp_ind)){
                len.of.loop <- 1:nrow(private$org_file)
            } else {
                len.of.loop <- self$which_exp_ind
            }

            # if k is not specified, perform calculation over c(2,4,6,8).
            if(length(k) == 1){
                private$k <- k
            } else {
                private$k <- c(4, 6, 8)
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
                private$breakpoint_path <- paste0("../../data/", 
                                                  private$bp_exp,
                                                  "/kmertone")
                private$ref_path <- paste0("../../data/ref/", 
                                           private$org_file[i, `Reference genome folder`])
                private$output_path <- paste0("../data/kmertone/", 
                                              private$bp_exp)

                for(kmer in private$k){
                    private$run_kmertone(kmer = kmer)
                    private$get_kmer_weight_factor(kmer = kmer)
                }
            }
        },

        #' @description
        #' Get heatmap tracts from kmertone data.
        #' @param k Numeric vector of kmer to perform calculation on.
        #' @return None.
        get_kmer_tracts = function(k){
            if(is.null(self$which_exp_ind)){
                len.of.loop <- 1:nrow(private$org_file)
            } else {
                len.of.loop <- self$which_exp_ind
            }

            # if k is not specified, perform calculation over c(2,4,6,8).
            if(length(k) == 1){
                private$k <- k
            } else {
                private$k <- c(4, 6, 8)
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
                private$breakpoint_path <- paste0("../../data/", 
                                                private$bp_exp,
                                                "/kmertone")
                private$ref_path <- paste0("../../data/ref/", 
                                        private$org_file[i, `Reference genome folder`])
                private$output_path <- paste0("../data/kmertone/", 
                                            private$bp_exp)

                for(kmer in private$k){
                    private$generate_table(kmer = kmer)
                    private$get_kmer_heatmaps(kmer = kmer)
                }
            }
        }
    ), 
    private = list(
        #' @field org_file Data.table of the full org_file.csv
        org_file = NULL,

        #' @field bp_exp Character vector of fragment type and experiment.
        bp_exp = NULL,

        #' @field breakpoint_path Character vector of path to kmertone folder.
        breakpoint_path = NULL,

        #' @field ref_path Character vector of path to reference genome folder.
        ref_path = NULL,

        #' @field output_path Character vector of path to save kmertone results.
        output_path = NULL,

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
        #' @param kmer Numeric vector of kmer to perform calculation on.
        #' @return None.
        generate_table = function(kmer){
            k.mers <- do.call(data.table::CJ, 
                              rep(list(c("A", "C", "G", "T")), kmer))
            self$kmer_list <- k.mers[, do.call(paste0, .SD)]
            
            rev.comp <- as.character(
                Biostrings::reverseComplement(Biostrings::DNAStringSet(self$kmer_list))
            )
            kmer_ref <- data.table('kmer' = self$kmer_list, 'rev.comp' = rev.comp)
            kmer_ref[, cond := ifelse(seq(1:nrow(.SD)) < match(kmer, rev.comp), 
            TRUE, ifelse(kmer == rev.comp, TRUE, FALSE))]
            self$kmer_ref <- kmer_ref[cond == TRUE, .(kmer)]
        },

        #' @description
        #' Run kmertone software
        #' @param kmer Numeric vector of kmer to perform calculation on.
        #' @return None.
        run_kmertone = function(kmer){
            # check for strand sensitivity
            sample.files <- list.files(
                path = private$breakpoint_path, 
                pattern = "csv|txt",
                full.names = TRUE
            )
            if(length(sample.files) > 0){
                dt <- fread(sample.files[1], nrows=1)
            }
            strand.sensitive <- ifelse("strand" %in% colnames(dt), TRUE, FALSE)

            # obtain long-range sequence context for control region cut-off
            dt.rmsd <- fread(paste0("../../03_FitCurves/data/ranges/", 
                                    "kmer_", kmer, "_Ranges_cutoffs_from_", 
                                    "clustering_all-exp.csv"))
            dt.rmsd <- dt.rmsd[exp == private$bp_exp, ]
            to.keep <- colnames(dt.rmsd)[grepl(pattern = "range", x = colnames(dt.rmsd))]
            longest.range <- ceiling(max(dt.rmsd[, ..to.keep], na.rm = TRUE)/2)
            print(paste0("Control region range = ", c(longest.range, longest.range+400)), 
                        quote = FALSE)

            # run kmertone
            dir.create(
                path = private$output_path,
                showWarnings = FALSE,
                recursive = TRUE
            )
            kmertone(
                case.coor.path=private$breakpoint_path, 
                genome.name="unknown", 
                strand.sensitive=strand.sensitive,
                k=kmer,
                ctrl.rel.pos=c(longest.range, longest.range+400), 
                case.pattern=NULL,
                output.path=private$output_path, 
                case.coor=NULL, 
                genome=NULL,
                genome.path=private$ref_path, 
                rm.case.kmer.overlaps=FALSE,
                case.length=2,
                merge.replicate=TRUE, 
                kmer.table=NULL,
                module="score", 
                ncpu=self$cores
            )
        },

        #' @description
        #' Get weight factors for kmers from kmertone
        #' @param kmer Numeric vector of kmer to perform calculation on.
        #' @return None.
        get_kmer_weight_factor = function(kmer){
            # load files
            files <- list.files(
                path = private$output_path, 
                pattern = paste0("score_", kmer)
            )
            files <- stringr::str_sort(files, numeric = TRUE)

            # import files
            data.set <- fread(
                file = files, 
                header = TRUE,
                showProgress = FALSE
            )

            # obtain weighting factors for use in breakpoint model
            norm.case <- data.set$case/sum(data.set$case, na.rm = TRUE)
            norm.control <- data.set$control/sum(data.set$control, na.rm = TRUE)
            ratio <- norm.case/norm.control 

            data.set <- data.table(
                kmer = data.set$kmer,
                weight.factor = ratio
            )

            dir.create(
                path = paste0("../data/weight_factor/", private$bp_exp),
                showWarnings = FALSE,
                recursive = TRUE
            )
            fwrite(
                data.set,
                file = paste0("../data/weight_factor/", 
                              private$bp_exp, 
                              "/kmer_", kmer, ".csv"), 
                row.names = FALSE
            )
        },

        #' @description
        #' Get heatmap tracts from kmertone data.
        #' @param kmer Numeric vector of kmer to perform calculation on.
        #' @return None.
        get_kmer_heatmaps = function(kmer){
            # import kmertone data
            all.files <- list.files(
                path = "../data/kmertone/", 
                pattern = paste0("score_", kmer, " *"), 
                recursive = TRUE,
                full.names = TRUE
            )
            all.files <- all.files[grepl(pattern = "true", x = all.files)]

            
        }
    )
)