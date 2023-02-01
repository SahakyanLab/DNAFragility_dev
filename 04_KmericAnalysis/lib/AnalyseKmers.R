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

        #' @field statistic Character vector of summation technique
        #' for a group of k-mer scores per breakage type.
        statistic = "mean",

        initialize = function(which_exp_ind, cores, statistic){
            if(!missing(which_exp_ind)) self$which_exp_ind <- which_exp_ind
            if(!missing(cores)) self$cores <- cores
            if(!missing(statistic)) self$statistic <- statistic

            # get full org_file.csv
            private$get_org_file()
        },

        #' @description
        #' Run kmertone software and obtain weight factors for each kmer.
        #' @param k Numeric vector of kmer to perform calculation on.
        #' @return None.
        get_kmer_enrichment = function(k = NULL){
            if(is.null(self$which_exp_ind)){
                self$which_exp_ind <- which(
                    private$org_file$`DSB Map` == "TRUE"
                )
            }
            len.of.loop <- self$which_exp_ind
            categories <- unique(private$org_file$Category_Main)
            private$categories <- sort(categories[which(categories != "")])            

            # if k is not specified, perform calculation over c(2,4,6,8).
            if(length(k) == 1){
                private$k <- k
            } else {
                private$k <- c(4, 6, 8)
            }

            # loop over each experiment
            for(i in len.of.loop){
                start.time <- Sys.time()

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
                private$rmsd_LR_threshold <- private$org_file[i, `RMSD?`]

                cur.msg <- paste0("Processing exp ", i, "/", 
                                  nrow(private$org_file),
                                  ". ", private$bp_exp)
                cat(cur.msg, "\n", sep = "")                                             

                for(kmer in private$k){
                    private$kmer <- kmer
                    private$run_kmertone()
                    private$get_kmer_weight_factor()
                }
            }
        },

        #' @description
        #' Get heatmap tracts from kmertone data.
        #' @param k Numeric vector of kmer to perform calculation on.
        #' @param action Character vector of "zscore" or "ratios".
        #' @param get_table Boolean. If TRUE, generates k-mer tables even if it exists.
        #' @return None.
        get_kmer_tracts = function(k = NULL, action = "zscore", get_table = FALSE){
            # if k is not specified, perform calculation over c(2,4,6,8).
            if(length(k) == 1){
                private$k <- k
            } else {
                private$k <- c(4, 6, 8)
            }

            if(dir.exists("../data/kmertone/QueryTable")){
                for(kmer in private$k){
                    private$kmer <- kmer
                    check.files <- list.files(
                        path = "../data/kmertone/QueryTable",
                        pattern = paste0("QueryTable_kmer-", 
                                         private$kmer, "_", action)
                    )
                    if(length(check.files) == 0 | get_table){
                        private$generate_table()
                        private$generate_querytable(action = action)
                    }
                    private$get_kmer_heatmaps(action = action)
                } 
            }
        }
    ), 
    private = list(
        #' @field org_file Data.table of the full org_file.csv
        org_file = NULL,

        #' @field categories Character vector of all breakage types.
        categories = NULL,

        #' @field bp_exp Character vector of fragment type and experiment.
        bp_exp = NULL,

        #' @field breakpoint_path Character vector of path to kmertone folder.
        breakpoint_path = NULL,

        #' @field ref_path Character vector of path to reference genome folder.
        ref_path = NULL,

        #' @field output_path Character vector of path to save kmertone results.
        output_path = NULL,

        #' @field rmsd_LR_threshold Boolean. If TRUE, will use long-range sequence
        #'  context influence for the control region of kmertone.
        rmsd_LR_threshold = FALSE,

        #' @field k Numeric vector of k-mer size for all iterations.
        k = NULL,

        #' @field kmer Numeric vector of k-mer size for current iteration.
        kmer = NULL,

        #' @description
        #' Import full org_file.csv and filter for rows to be processed.
        #' @return None.
        get_org_file = function(){
            private$org_file <- fread(
                "../../data/org_file.csv",
                showProgress = FALSE
            )
            private$org_file[, exp := paste0(`Fragmentation type`, "/", 
                                             `Experiment folder`)]                        
        },

        #' @description 
        #' A utility function to generate k-mers.
        #' Then, only keeps first occurring k-mer in lexicological order.
        #' @return None.
        generate_table = function(){
            k.mers <- do.call(data.table::CJ, 
                              rep(list(c("A", "C", "G", "T")), private$kmer))
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
        #' Generates k-mer probability ratios or z-score query table.
        #' @return None.
        generate_querytable = function(action){
            # progress message
            t1 <- Sys.time()
            cur.msg <- paste0("Generating ", action, " ", 
                              private$kmer, "-mer table")
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # Query_Table for k-mer breakage propensities
            df <- lapply(private$categories, function(category){
                return(private$get_data(
                    category = category, 
                    action = action
                ))
            })
            df <- rbindlist(df)
            df <- as_tibble(df) %>% 
                tidyr::spread(key = "kmer", value = "value")

            dir.create(
                path = "../data/kmertone/QueryTable/",
                showWarnings = FALSE,
                recursive = TRUE
            )
            fwrite(
                x = df, 
                file = paste0("../data/kmertone/QueryTable/",
                              "QueryTable_kmer-", private$kmer, 
                              "_", action, ".csv")
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' A utility function to calculate summarised probability ratios 
        #' or z-scores per group of breakage source.
        #' @param category Character vector of a single breakage source.
        #' @return A data.table of k-mer scores per breakage category.
        get_data = function(category, action){
            org.filtered <- private$org_file[Category_Main == category, exp]

            # mean/median of same breakage type
            out <- lapply(1:length(org.filtered), function(x){
                return(private$summarise_score(
                    file = org.filtered[x],
                    action = action
                ))
            })

            out <- switch(self$statistic,
                "mean" = {
                    sapply(1:length(out[[1]]), function(x){
                        mean(sapply(out, `[[`, x), na.rm = TRUE)
                    })
                },
                "median" = {
                    sapply(1:length(out[[1]]), function(x){
                        median(sapply(out, `[[`, x), na.rm = TRUE)
                    })
                }
            )

            datatable <- data.table(
                kmer = self$kmer_ref$kmer,
                value = out,
                category = as.factor(category)
            )
            return(datatable)
        },

        #' @description 
        #' A utility function to calculate k-mer probability ratios or z-scores.
        #' @param file Character vector of the specified breakage type.
        #' @return A data.table of scores.
        summarise_score = function(file, action){
            files <- list.files(
                path = paste0("../data/kmertone/", file),
                pattern = paste0("score_", private$kmer),
                full.names = TRUE
            )

            data.set <- fread(
                file = files,
                showProgress = FALSE,
                select = c("case", "control", "z")
            )

            data.set[, kmer := self$kmer_list]            
            data.set <- data.set[!which(is.na(data.set$z)), ]
            data.set <- data.set[!which(is.infinite(data.set$z)), ]

            # only keep first occurring kmer in lexicological order 
            data.set <- data.set[match(self$kmer_ref$kmer, data.set$kmer)]            
            output <- switch(action,
                "ratio" = {
                    norm.case <- data.set$case/sum(data.set$case, 
                                                   na.rm = TRUE)
                    norm.control <- data.set$control/sum(data.set$control, 
                                                         na.rm = TRUE)
                    ratio <- norm.case/norm.control 
                },
                "zscore" = data.set$z
            )
            return(output)
        },

        #' @description
        #' Run kmertone software
        #' @return None.
        run_kmertone = function(){
            # check for strand sensitivity
            sample.files <- list.files(
                path = private$breakpoint_path, 
                pattern = ".csv",
                full.names = TRUE
            )
            if(length(sample.files) > 0){
                dt <- fread(sample.files[1], nrows=1)
            }
            strand.sensitive <- ifelse("strand" %in% colnames(dt), TRUE, FALSE)

            if(private$rmsd_LR_threshold){
                # obtain long-range sequence context for control region cut-off
                dt.rmsd <- fread(paste0("../../03_FitCurves/data/ranges/", 
                                        "kmer_", private$kmer, "_Ranges_cutoffs_from_", 
                                        "clustering_all-exp.csv"))
                dt.rmsd <- dt.rmsd[exp == private$bp_exp, ]
                to.keep <- colnames(dt.rmsd)[grepl(pattern = "range", x = colnames(dt.rmsd))]
                longest.range <- ceiling(max(dt.rmsd[, ..to.keep], na.rm = TRUE)/2)
                control.regions <- c(longest.range, longest.range+1000)
                case.len <- 2
            } else {
                # take longest range
                files <- list.files(
                    path = paste0("../../data/", private$bp_exp),
                    pattern = "bed|narrowpeak|broadpeak|csv",
                    full.names = TRUE,
                    ignore.case = TRUE
                )
                pos <- regexpr("\\.([[:alnum:]]+)$", files)
                file.ext <- ifelse(pos > -1, 
                    substring(text = files, first = pos+1), 
                    ""
                )

                if(grepl(pattern = "csv", x = file.ext, ignore.case = TRUE)){
                    df <- fread(files, showProgress = FALSE)
                    setnames(df, c("V1", "V2", "V3"))
                    df <- df[, .(V1, V2, V3)]
                    setnames(df, c("seqnames", "start", "end"))
                    df <- plyranges::as_granges(df)            
                } else if(grepl(pattern = "narrowpeak|bed|broadpeak", x = file.ext, ignore.case = TRUE)){
                    if(grepl(pattern = "narrowpeak", x = file.ext, ignore.case = TRUE)){
                        df <- plyranges::read_narrowpeaks(files)
                    } else if(grepl(pattern = "bed|broadpeak", x = file.ext, ignore.case = TRUE)){
                        df <- tryCatch({
                            plyranges::read_bed(files)
                        }, error = function(e){
                            temp <- fread(files, showProgress = FALSE)
                            temp <- temp[, .(V1, V2, V3)]
                            setnames(temp, c("seqnames", "start", "end"))
                            temp[, start := start+1] # bed-files are 0-indexed
                            return(plyranges::as_granges(temp))
                        })
                    }
                    if(ncol(attr(df, "elementMetadata")) > 0){
                        df <- df %>% dplyr::select(-colnames(attr(df, "elementMetadata")))
                    }
                }

                df <- df %>% dplyr::filter(seqnames == "chr1")
                # control start pos: average of top 5% of widths
                top.five.perc <- sort(width(ranges(df)), decreasing = TRUE)
                top.five.perc <- top.five.perc[1:ceiling(length(top.five.perc)*0.05)]
                top.five.perc <- ceiling(mean(top.five.perc, na.rm = TRUE))
                # control end pos: 2x of max width or 1000 bp for consistency with RMSD cutoffs
                upper.control <- max(c(ceiling(max(width(ranges(df)), na.rm = TRUE)*2), 1000))
                control.regions <- c(top.five.perc, upper.control)
                case.len <- NULL
            }                        
            cur.msg <- paste0("Control region range: [", control.regions[1], 
                              ",", control.regions[2], "]")
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")                             

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
                k=private$kmer,
                ctrl.rel.pos=control.regions, 
                case.pattern=NULL,
                output.path=private$output_path, 
                case.coor=NULL, 
                genome=NULL,
                genome.path=private$ref_path, 
                rm.case.kmer.overlaps=FALSE,
                case.length=case.len,
                merge.replicate=TRUE, 
                kmer.table=NULL,
                module="score", 
                ncpu=self$cores
            )
        },

        #' @description
        #' Get weight factors for kmers from kmertone
        #' @return None.
        get_kmer_weight_factor = function(){
            # load files
            files <- list.files(
                path = private$output_path, 
                pattern = paste0("score_", private$kmer),
		        full.names = TRUE
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
                              private$bp_exp, "/kmer_", 
                              private$kmer, ".csv"), 
                row.names = FALSE
            )
        },

        #' @description
        #' Get heatmap tracts from kmertone data.
        #' @param action Character vector of "zscore" or "ratios".
        #' @return None.
        get_kmer_heatmaps = function(action){
            t1 <- Sys.time()
            cur.msg <- paste0("Clustering heatmaps for kmer: ", private$kmer)
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # import processed kmertone data
            data.set <- fread(
                file = paste0("../data/kmertone/QueryTable/QueryTable_kmer-", 
                              private$kmer, "_", action, ".csv"),
                showProgress = FALSE
            )
            data.set.labels <- data.set$category
            data.set.kmers <- colnames(data.set[,-"category"])
            
            # normalise by row (experiment)
            data.set <- apply(data.set[,-"category"], 1, scale, 
                              scale = TRUE, center = TRUE)
            data.set <- t(data.set)
            rownames(data.set) <- data.set.labels
            colnames(data.set) <- data.set.kmers

            # hierarchical clustering
            rownames(data.set) <- data.set.labels
            df.dist <- dist(1-data.set) %>% suppressWarnings()
            df.dendro <- as.dendrogram(hclust(df.dist, method = "ward.D"))
            transpose.df <- t(data.set)
            df.dendro.transpose <- as.dendrogram(hclust(dist(1-transpose.df), method = "ward.D"))
            cols <- 10
            CI <- mean(data.set, na.rm = TRUE)+c(-1, 1)*1.96*sd(data.set, na.rm = TRUE)

            dir.create(
                path = "../figures/kmertone_tracts",
                showWarnings = FALSE,
                recursive = TRUE
            )
            pdf(
                paste0("../figures/kmertone_tracts/", private$kmer, 
                       "-mer_", action, "-clustering.pdf"), 
                height = 7, 
                width = 12
            )
            heatmap.2(
                data.set,
                offsetRow = 0,
                offsetCol = 0,
                Rowv = df.dendro,
                Colv = df.dendro.transpose,
                dendrogram = "both",
                revC = FALSE,
                trace = "none",
                density.info = "none",
                col = hcl.colors(cols, "RdYlGn"),
                breaks = seq(CI[1], CI[2], length.out = cols+1),
                notecol = "black",
                cexRow = 0.7,
                labRow = rownames(data.set),
                labCol = "",
                margins = c(3,19),
                symkey = TRUE,
                key.xlab = "Scaled z-scores"
            )
            plot.save <- dev.off()

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        }
    )
)
