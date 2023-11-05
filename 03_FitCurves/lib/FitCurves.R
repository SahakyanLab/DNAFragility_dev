FitCurves <- R6::R6Class(
    classname = "FitCurves",
    public = list(
        #' @field chr Numeric vector of chromosome number.
        chr = NULL,

        #' @field which_exp_ind Numeric vector of the index of the experiment 
        #' to process from org_file. If NULL, all exp will be processed.
        which_exp_ind = NULL,

        #' @field results List of results stored if return_vals is TRUE.
        results = NULL,

        initialize = function(chr, which_exp_ind, from_cluster, return_vals, rmsd_values){
            if(!missing(chr)) self$chr <- chr
            if(!missing(which_exp_ind)) self$which_exp_ind <- which_exp_ind
            if(!missing(from_cluster)) private$from_cluster <- from_cluster
            if(!missing(return_vals)) private$return_vals <- return_vals
            
            # get full org_file.csv
            private$get_org_file()

            self$results <- vector(
                mode = "list", 
                length = length(self$chr)
            )
            names(self$results) <- paste0("chr", self$chr)

            if(!missing(rmsd_values)){
                for(i in 1:length(rmsd_values)){
                    private$rmsd_values[[names(rmsd_values)[i]]] <- rmsd_values[[names(rmsd_values)[i]]]
                }
            }

            if(is.null(private$rmsd_values)){
                if(is.null(self$which_exp_ind)){
                    self$which_exp_ind <- which(
                        (private$org_file$`DSB Map` == "TRUE") & 
                        (private$org_file$`RMSD?` == "TRUE")
                    )
                    # if(!private$from_cluster) private$to_cluster_curves <- TRUE
                }
                private$len_of_loop <- self$which_exp_ind
            } 
        },

        #' @description
        #' Generate plots from RMSD calculations for each experiment.
        #' @param k size of kmer to plot and fit Gaussian curves to.
        #' @param annot_plots if TRUE, will add 95% C.I. of each curve onto plot.
        #' @return None.
        generate_rmsd_plots = function(k = c(2,4,6,8), annot_plots = TRUE){
            start.time <- Sys.time()
            cur.msg <- "Generating RMSD plots and fitting Gaussian curves"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")

            if(private$from_cluster) plot.vec <- vector(mode = "list")

            if(is.null(private$rmsd_values)){
                for(i in private$len_of_loop){
                    private$bp_exp <- paste0(
                        private$org_file[i, `Fragmentation type`], "/",
                        private$org_file[i, `Experiment folder`]
                    )
                    private$category <- private$org_file[i, Category_Main]
                    private$category_general <- private$org_file[i, Category_general]

                    cur.msg <- paste0("Processing ", i, "/", 
                                    nrow(private$org_file),
                                    ". ", private$bp_exp)
                    cat(cur.msg, "\n", sep = "")

                    # if k is not specified, perform calculations on c(2,4,6,8)
                    if(length(k) == 1){
                        if(k == 8) private$fits <- c("kmer_8" = 3)
                    } else {
                        private$fits <- c(
                            "kmer_2" = 3, 
                            "kmer_4" = 3, 
                            "kmer_6" = 3, 
                            "kmer_8" = 3
                        )
                    }

                    for(chr in self$chr){
                        t1 <- Sys.time()
                        cur.msg <- paste0("Calculating RMSD values for chromosome ", chr)
                        l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                        cat(cur.msg, l, sep = "")

                        if(!private$from_cluster){
                            private$plot_rmsd(
                                annot_plots = annot_plots,
                                chr = chr
                            )

                            total.time <- Sys.time() - t1
                            cat("DONE! --", signif(total.time[[1]], 2), 
                                attr(total.time, "units"), "\n")
                            next
                        }

                        # if data from averaged summed sequence context dependence plots,
                        # save the combined plots
                        plot.vec <- c(
                            plot.vec, 
                            private$plot_rmsd(annot_plots = annot_plots)
                        )

                        total.time <- Sys.time() - t1
                        cat("DONE! --", signif(total.time[[1]], 2), 
                            attr(total.time, "units"), "\n")
                    }
                }
            } else {
                for(chr in self$chr){
                    t1 <- Sys.time()
                    cur.msg <- paste0("Calculating RMSD values for chromosome ", chr)
                    l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                    cat(cur.msg, l, sep = "")

                    if(k == 8){
                        private$fits <- c("kmer_8" = 3)
                    } else {
                        stop("NotImplementedError.")
                    }
                    private$plot_rmsd(
                        annot_plots = annot_plots,
                        chr = chr
                    )

                    total.time <- Sys.time() - t1
                    cat("DONE! --", signif(total.time[[1]], 2), 
                        attr(total.time, "units"), "\n")
                }
            }
            
            if(private$from_cluster){
                pdf(
                    file = paste0("../figures/", private$org_file[i, `Fragmentation type`], 
                                "/chr", self$chr, "_RMSD_all_GMMFit.pdf"),
                    height = 20, width = 7
                )
                p1 <- do.call(gridExtra::grid.arrange, c(plot.vec, ncol = 1))
                pic.saved <- dev.off()                
            }

            if(private$to_cluster_curves){
                # cluster rmsd tracts into ranges
                for(kmer in c(2,4,6,8)){                   
                    private$cluster_curves_into_ranges(kmer)
                }
            }

            # if(is.null(self$which_exp_ind)){
            #     # clean up the last 3 columns if not empty
            #     system("/usr/local/bin/bash check_to_clean.sh")

            #     # init heatmap tracts from RMSD data
            #     for(kmer in c(4,6,8)){
            #         private$get_rmsd_tracts(kmer = kmer)
            #     }

            #     # hierarchical clustering to find optimal curve fits
            #     for(kmer in c(4,6,8)){
            #         private$save_curve_counts(kmer = kmer)
            #     }

            #     private$auto_fit <- FALSE
            #     select.cols <- c("kmer_4", "kmer_6", "kmer_8")
            #     private$fits <- setNames(
            #         as.numeric(private$org_file[i, ..select.cols]),
            #         colnames(private$org_file[i, ..select.cols])
            #     )
            #     private$plot_rmsd()

            #     # update heatmap tracts from RMSD data
            #     for(kmer in c(4,6,8)){
            #         private$get_rmsd_tracts(kmer = kmer)
            #     }
            # }

            # time taken for full processing for this experiment
            final.t <- Sys.time() - start.time
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cat("Final time taken:", signif(final.t[[1]], digits = 2), 
                attr(final.t, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
        }
    ), 
    private = list(
        #' @field org_file Data.table of the full org_file.csv
        org_file = NULL,

        #' @field bp_exp Character vector of fragment type and experiment.
        bp_exp = NULL,

        #' @field category Character vector of the breakage category.
        category = NULL,

        #' @field category_general Character vector of the breakage sub-category.
        category_general = NULL,

        #' @field auto_fit Boolean. If TRUE, will not optimise 
        #' number of curves to fit.
        auto_fit = TRUE,

        #' @field len_of_loop Numeric vector of how many experiments to process.
        len_of_loop = NULL,

        #' @field fits Named vector of curves to fit per k-mer RMSD values.
        fits = NULL,

        #' @field to_cluster_curves Boolean. If TRUE, will cluster curves into ranges.
        to_cluster_curves = FALSE,

        #' @field from_cluster Boolean. If TRUE, will generat plots and fit curves
        #' from averaged summed sequence context depedence datasets.
        from_cluster = FALSE,

        #' @field return_vals Boolean. If TRUE, will return values instead of saving as file.
        return_vals = FALSE,

        #' @field rmsd_values Numeric vector of RMSD values. 
        rmsd_values = NULL,

        #' @description
        #' Import full org_file.csv and filter for rows to be processed.
        #' @return None.
        get_org_file = function(){
            if(private$from_cluster){
                private$org_file <- fread(
                    "../../data/org_file_from_clustering.csv",
                    showProgress = FALSE
                )                
            } else {
                private$org_file <- fread(
                    "../../data/org_file.csv",
                    showProgress = FALSE
                )
            }
        },

        #' @description
        #' Generate RMSD plots for all k-mers available
        #' @param annot_plots if TRUE, will add 95% C.I. of each curve onto plot.
        #' @param chr Numeric vector of chromosome number.
        #' @return None.
        plot_rmsd = function(annot_plots, chr){
            if(is.null(private$rmsd_values)){
                # load data sets
                files <- list.files(
                    path = paste0("../data/", private$bp_exp),
                    pattern = "^rmsd_kmer_",
                    full.names = TRUE
                )
                files <- stringr::str_sort(files, numeric = TRUE)
                file.name <- stringr::str_extract(
                    string = basename(files),
                    pattern = "kmer_[:digit:]"
                )

                data.sets <- lapply(1:length(files), function(x){
                    out <- readRDS(file = files[x])
                    limits <- length(out)/2-1
                    out <- as_tibble(out) %>% 
                        dplyr::mutate(
                            kmer = as.factor(file.name[x]), 
                            x = -limits:(length(out)-limits-1)) %>% 
                        dplyr::rename(y = value)
                    return(out)
                })
                data.sets <- do.call(rbind, data.sets)

                # Overall RMSD plots
                plots <- data.sets %>%
                    dplyr::mutate(kmer = 
                        stringr::str_replace_all(
                            string = stringr::str_to_title(kmer),
                            pattern = "_",
                            replacement = " "
                        )) %>%              
                    ggplot(aes(x = x, y = y)) + 
                        geom_line(linewidth = 0.8) + 
                        facet_wrap(~kmer, ncol = 4, scales = "free_y") + 
                        theme_bw() + 
                        coord_cartesian(ylim = c(0, NA)) + 
                        theme(text = element_text(size = 25)) + 
                        labs(
                            x = "Position away from breakpoint",
                            y = "RMSD"
                        )

                dir.create(
                    path = paste0("../figures/", private$bp_exp),
                    showWarnings = FALSE,
                    recursive = TRUE
                )

                height <- ifelse(length(files) == 1, 7, 6)
                width <- ifelse(length(files) == 1, 7, 21)
                if(!private$return_vals){
                    ggsave(
                        filename = paste0("../figures/", private$bp_exp, 
                                            "/chr", chr, 
                                            "_sidebyside_RMSD.pdf"),
                        plot = plots, 
                        height = height,
                        width = width
                    )
                }
            }

            # GMM plots
            if(is.null(private$rmsd_values)){
                data.sets <- data.sets %>%
                    dplyr::mutate(dplyr::across(where(is.factor), as.character)) %>% 
                    dplyr::filter(
                        kmer == "kmer_2" |
                        kmer == "kmer_4" |
                        kmer == "kmer_6" |
                        kmer == "kmer_8"
                    )
            } else {
                limits <- length(private$rmsd_values[[paste0("chr", chr)]])/2-1
                data.sets <- tibble(
                    x = -limits:(length(private$rmsd_values[[paste0("chr", chr)]])-limits-1),
                    y = private$rmsd_values[[paste0("chr", chr)]],
                    kmer = "kmer_8"
                )
            }

            p <- vector(mode = "list", length = length(private$fits))
            for(k in 1:length(private$fits)){
                if(!is.na(private$fits[k])){
                    dat <- data.sets %>% 
                        dplyr::filter(kmer == names(private$fits[k]))

                    curvefits <- private$fit_gmm(
                        dat = dat, 
                        ind = names(private$fits[k]), 
                        nr.of.curves = private$fits[k]
                    )

                    if(class(curvefits[length(curvefits)][[1]]) != "nls"){
                        private$fits[k] <- 2
                        curvefits <- private$fit_gmm(
                            dat = dat, 
                            ind = names(private$fits[k]), 
                            nr.of.curves = private$fits[k]
                        )   
                    }

                    if(class(curvefits[length(curvefits)][[1]]) != "nls"){
                        private$fits[k] <- 1
                        curvefits <- private$fit_gmm(
                            dat = dat, 
                            ind = names(private$fits[k]), 
                            nr.of.curves = private$fits[k]
                        )   
                    }

                    output <- private$make_plot(
                        dat = dat, 
                        k = names(private$fits[k]), 
                        curve.vals = curvefits, 
                        nr.of.curves = private$fits[k],
                        annot_plots = annot_plots
                    )
                    p[[k]] <- output[[1]]

                    # percent contribution of each gaussian curve towards breakage
                    df <- as_tibble(curvefits[[length(curvefits)-3]])*100
                    if(ncol(df) == 1) df <- dplyr::mutate(df, curve.two = NA_real_)
                    if(ncol(df) == 2) df <- dplyr::mutate(df, curve.three = NA_real_)

                    # range of influence based on 95 percent confidence intervals
                    output.CIlst <- output[[2]]
                    curves <- output.CIlst[seq(2, length(output.CIlst), 2)]-
                              output.CIlst[seq(1, length(output.CIlst), 2)]
                    curves <- c(curves, rep(NA_real_, 3-length(curves)))

                    # sd of curves
                    SD <- as.numeric(curvefits[[length(curvefits)-2]])
                    SD <- c(SD, rep(NA_real_, 3-length(SD)))

                    # contributino of the peak intensity of each gaussian curves
                    peak.intensity <- as.numeric(curvefits[[length(curvefits)-1]])*100
                    peak.intensity <- c(peak.intensity, rep(NA_real_, 3-length(peak.intensity)))

                    # combine results
                    df <- rbind(df, curves, SD, peak.intensity)
                    df <- df %>% 
                        dplyr::mutate(
                            exp = private$bp_exp,
                            category = private$category,
                            category_general = private$category_general,
                            rowid = c("contribution", "ranges", "SD", "peak.intensity"), 
                            .before = 1,
                        )

                    # save the individual points of the gaussian curve
                    curve.vals <- lapply(1:unname(private$fits[k]), function(x){
                        temp <- data.table(vals = curvefits[[x]])  

                        if(x == 1){
                            col.label <- "curve.one"
                        } else if(x == 2){
                            col.label <- "curve.two"
                        } else if(x == 3){
                            col.label <- "curve.three"
                        }
                        
                        temp[, curve := col.label]
                        return(temp)
                    })
                    curve.vals <- rbindlist(curve.vals)

                    if(!private$auto_fit){
                        # extract cut-off values for 
                        cutoff.df <- fread("../data/Ranges_cutoffs.csv")
                        cutoff.df <- as_tibble(cutoff.df) %>% 
                            dplyr::pull(Cluster, names(private$fits[k]))

                        temp <- df %>% 
                            dplyr::filter(rowid == "ranges") %>% 
                            dplyr::select(contains("curve")) %>% 
                            tidyr::gather(key = "ranges", value = "curves") %>% 
                            dplyr::pull(ranges, curves)

                        # short range
                        short.range <- which(
                            as.integer(names(temp)) <= as.integer(names(cutoff.df[3]))
                        )
                        if(any(short.range)){
                            if(length(short.range) > 1){
                            # if multiple curves fall within the same range
                            # then at least 1 is a redundant curve to be removed
                            # after this, update the curve number on org_file.csv
                            contribution.temp <- df %>% 
                                dplyr::filter(rowid == "contribution") %>% 
                                dplyr::select(-c(1:3)) %>% 
                                tidyr::gather(key = "curves", value = "contribution") %>% 
                                dplyr::pull(contribution, curves)

                            colnames(df)[4:6][short.range] <- "short.range"
                            } else {
                            temp[short.range] <- cutoff.df[3]
                            colnames(df)[4:6][short.range] <- cutoff.df[3]
                            }
                        }

                        # mid range
                        mid.range <- which(
                            (as.integer(names(temp)) > as.integer(names(cutoff.df[3]))) & 
                            (as.integer(names(temp)) <= as.integer(names(cutoff.df[2])))
                        )
                        if(any(mid.range)){
                            if(length(mid.range) > 1){
                                # if multiple curves fall within the same range
                                # then at least 1 is a redundant curve to be removed
                                # after this, update the curve number on org_file.csv
                                contribution.temp <- df %>% 
                                    dplyr::filter(rowid == "contribution") %>% 
                                    dplyr::select(-c(1:3)) %>% 
                                    tidyr::gather(key = "curves", value = "contribution") %>% 
                                    dplyr::pull(contribution, curves)

                                colnames(df)[4:6][mid.range] <- "mid.range"
                            } else {
                                temp[mid.range] <- cutoff.df[2]
                                colnames(df)[4:6][mid.range] <- cutoff.df[2]
                            }
                        }

                        # long range
                        long.range <- which(
                            as.integer(names(temp)) > as.integer(names(cutoff.df[2]))
                        )
                        if(any(long.range)){
                            if(length(long.range) > 1){
                                # if multiple curves fall within the same range
                                # then at least 1 is a redundant curve to be removed
                                # after this, update the curve number on org_file.csv
                                contribution.temp <- df %>% 
                                    dplyr::filter(rowid == "contribution") %>% 
                                    dplyr::select(-c(1:3)) %>% 
                                    tidyr::gather(key = "curves", value = "contribution") %>% 
                                    dplyr::pull(contribution, curves)

                                colnames(df)[4:6][long.range] <- "long.range"
                            } else {
                                temp[long.range] <- "long.range"
                                colnames(df)[4:6][long.range] <- "long.range"
                            }
                        }
                    }                    

                    if(is.null(private$rmsd_values)){
                        df <- df[, 1:(4+unname(private$fits[k]))]
                    }
                    all.gc <- as.data.table(output[[3]])

                    if(private$return_vals){
                        # gaussian curve stats
                        self$results[[paste0("chr", chr)]]$df <- df 

                        # curve values
                        self$results[[paste0("chr", chr)]]$curve_vals <- curve.vals

                        # linear combination of all gaussian curves
                        self$results[[paste0("chr", chr)]]$all_gc <- all.gc
                    } else {
                        dir.create(
                            path = paste0("../data/", private$bp_exp),
                            showWarnings = FALSE,
                            recursive = TRUE
                        )

                        # save the gaussian curves statistics
                        fwrite(
                            df,
                            file = paste0("../data/", private$bp_exp, "/key",
                                        "_stats_", names(private$fits[k]), ".csv"),
                            row.names = FALSE
                        )

                        # save the actual gaussian curves values
                        fwrite(
                            curve.vals,
                            file = paste0("../data/", private$bp_exp, 
                                        "/gaussian_curve_vals_",
                                        names(private$fits[k]), ".csv"),
                            row.names = FALSE
                        )

                        # save the linear combination of all gaussian curves values
                        fwrite(
                            all.gc,
                            file = paste0("../data/", private$bp_exp, 
                                        "/gaussian_curve_vals_linear_combinations_",
                                        names(private$fits[k]), ".csv"),
                            row.names = FALSE
                        )                    
                    }
                }
                # if(k == 1) print(p) # temp. avoid an unsolved ggplot bug
            }

            if(is.null(private$rmsd_values)){
                height <- ifelse(length(files) == 1, 5, 12)
                width <- ifelse(length(files) == 1, 7, 5)
                if(!private$return_vals){
                    pdf(
                        file = paste0("../figures/", private$bp_exp, 
                                    "/chr", chr, "_RMSD_all_GMMFit.pdf"),
                        height = height, width = width
                    )
                    p1 <- do.call(gridExtra::grid.arrange, c(p, ncol = 1))
                    pic.saved <- dev.off()
                }
            }

            if(private$from_cluster) return(p)
        },

        #' @description
        #' Fit gaussian curves to underlying RMSD values.
        #' @param dat Tibble of the RMSD values per k-mer.
        #' @param ind Numeric vector of the index.
        #' @param C.value Numeric vector of the starting gaussian coefficient.
        #' @param nr.of.curves Numeric vector of the number of curves to fit.
        #' @return None.
        fit_gmm = function(dat, ind, C.value = 10, nr.of.curves){
            # Fit to a model consisting of a pair of Gaussians. Note that nls is kind
            # of fussy about choosing "reasonable" starting guesses for the parameters.
            # It's even fussier if you use the default algorithm (i.e., Gauss-Newton)
            x <- dplyr::pull(dat, x)
            y <- dplyr::pull(dat, y)
            dat <- data.frame(x = x, y = y)

            # average normalised behaviour
            dat.norm <- tibble(
                x = x,
                y_fwd = scale(y, center = TRUE, scale = TRUE),
                y_rev = rev(y_fwd)) 
            dat.norm <- data.frame(
                x = x,
                y = rowMeans(dat.norm[2:3], na.rm = TRUE)
            )
            dat.norm <- dat.norm[1:(ceiling(nrow(dat.norm)/2)),]

            # grid-search based approach of starting values for nls convergence
            grid_search <- data.frame(
                sd1 = seq(1, 20, 0.5),
                sd2 = seq(4, 62, 1.5),
                sd3 = seq(10, 68, 1.5)                
            )
            grid_search$rss <- rep(NA, length(grid_search$sd1))
            
            nls.res <- lapply(1:nrow(grid_search), function(index){
                if(nr.of.curves == 1){
                    try_fit <- function(C.value = 10, sd1){
                        nls(
                        y ~ (C1 * exp(-(x-0)^2/(2 * sigma1^2)) +
                            min(y)),
                        data=dat,
                        start=list(C1=C.value, sigma1=sd1),
                        control = nls.control(
                            maxiter=50000, 
                            tol=1e-05, 
                            warnOnly = TRUE
                        ),
                        lower=rep(c(0, 0)),
                        upper=rep(c(NULL, 500)),
                        # upper=rep(c(NULL, Inf)),
                        trace=FALSE,
                        algorithm="port") %>% 
                        suppressWarnings()
                    }
                    fit <- tryCatch({
                        try_fit(
                            C.value = 10,
                            sd1 = grid_search$sd1[index]
                        )
                    }, error = function(e) return(2))
                } else if(nr.of.curves == 2){
                    try_fit <- function(C.value = 10, sd1, sd2){
                        nls(
                        y ~ (C1 * exp(-(x-0)^2/(2 * sigma1^2)) +
                                C2 * exp(-(x-0)^2/(2 * sigma2^2)) +
                            min(y)),
                        data=dat,
                        start=list(C1=C.value, sigma1=sd1,
                                    C2=C.value, sigma2=sd2
                                    ),
                        control = nls.control(
                            maxiter=50000, 
                            tol=1e-05, 
                            warnOnly = TRUE
                        ),
                        lower=rep(c(0, 0), 2),
                        # upper=rep(c(NULL, Inf), 3),
                        upper=rep(c(NULL, 500), 3),
                        trace=FALSE,
                        algorithm="port") %>% 
                        suppressWarnings()
                    }
                    fit <- tryCatch({
                        try_fit(
                            C.value = 10,
                            sd1 = grid_search$sd1[index],
                            sd2 = grid_search$sd2[index]
                        )
                    }, error = function(e) return(2))
                } else if(nr.of.curves == 3){
                    try_fit <- function(C.value = 10, sd1, sd2, sd3){
                        nls(
                        y ~ (C1 * exp(-(x-0)^2/(2 * sigma1^2)) +
                                C2 * exp(-(x-0)^2/(2 * sigma2^2)) +
                                C3 * exp(-(x-0)^2/(2 * sigma3^2)) +
                            min(y)),
                        data=dat,
                        start=list(C1=C.value, sigma1=sd1,
                                    C2=C.value, sigma2=sd2,
                                    C3=C.value, sigma3=sd3
                                    ),
                        control = nls.control(
                            maxiter=50000, 
                            tol=1e-05, 
                            warnOnly = TRUE
                        ),
                        lower=rep(c(0, 0), 3),
                        upper=rep(c(NULL, 500), 3),
                        # upper=rep(c(NULL, Inf), 3),
                        trace=FALSE,
                        algorithm="port") %>% 
                        suppressWarnings()
                    }
                    fit <- tryCatch({
                        try_fit(
                            C.value = 10,
                            sd1 = grid_search$sd1[index],
                            sd2 = grid_search$sd2[index],
                            sd3 = grid_search$sd3[index]
                        )
                    }, error = function(e) return(2))
                }

                if(class(fit) == "nls"){
                    check.params <- fit$m$getAllPars()
                    # check for zero coefficients
                    any.zero.coefs <- any(check.params[grepl(
                        pattern = "C", x = names(check.params)
                    )] == 0)
                    # check for any out-of-bound sigmas 
                    any.oob.sigmas <- any(check.params[grepl(
                        pattern = "sigma", x = names(check.params)
                    )] <= 1)
                    if(any.zero.coefs | any.oob.sigmas){
                        return(NULL)
                    } else {
                        grid_search$rss[index] <<- signif(
                            fit$m$deviance(), 
                            digits = 10
                        )
                        return(fit)
                    }
                }
                return(NULL)
            })

            # find best fit if exists
            if(all(is.na(grid_search$rss))){
                fit <- 2
            } else {
                best.fit.index <- which.min(grid_search$rss)
                fit <- nls.res[[best.fit.index]]
            }

            # error capture
            init.list <- vector(mode = "list", length = 4)
            init.list[[4]] <- fit
            if(class(fit) != "nls") return(init.list)

            # extract each parameter
            summary.fit.params <- fit$m$getAllPars()

            # extract rss
            rss <- sum(residuals(fit)^2)

            # # check for zero coefficients
            # any.zero.coefs <- any(summary.fit.params[grepl(
            #     pattern = "C", x = names(summary.fit.params)
            # )] == 0)
            # if(any.zero.coefs) return(init.list[[4]] <- 2)
            
            # fit each Gaussian curve separately to data
            draw.from.gaussian <- function(xs, C, SD){
                return(C*exp(-(xs-0)^2/(2*SD^2))+min(dat$y, na.rm = TRUE))
            }
            
            # fit curve 1
            fit_1 <- draw.from.gaussian(
                xs = x,
                C = summary.fit.params["C1"],
                SD = summary.fit.params["sigma1"]
            )
            peak.intensity.fit_1 <- max(fit_1, na.rm = TRUE)

            if(nr.of.curves == 1){
                return(
                    list(
                        fit_1, 
                        list("curve.one" = 1),
                        list("curve.one" = summary.fit.params["sigma1"]),
                        list("curve.one" = peak.intensity.fit_1),
                        # list("curve.one" = summary.fit.params["C1"]),                        
                        init.list[length(init.list)][[1]]
                    )
                )
            } else if(nr.of.curves > 1){
                # fit curve 2
                fit_2 <- draw.from.gaussian(
                    xs = x,
                    C = summary.fit.params["C2"],
                    SD = summary.fit.params["sigma2"]
                )
                peak.intensity.fit_2 <- max(fit_2, na.rm = TRUE)-min(dat$y, na.rm = TRUE)

                # calculate area under each curve
                integral.curve.one <- summary.fit.params["C1"]*summary.fit.params["sigma1"]*sqrt(2*pi)
                integral.curve.two <- summary.fit.params["C2"]*summary.fit.params["sigma2"]*sqrt(2*pi)

                if(nr.of.curves == 2){
                    # Calculate total area
                    integral.all.curves <- integral.curve.one+integral.curve.two

                    # percentage contribution of each curve
                    curve.one.contribution <- integral.curve.one/integral.all.curves
                    curve.two.contribution <- 1-curve.one.contribution

                    # peak intensity contribution of each gaussian curve normalised
                    peak.intensity.all <- (fit_1+fit_2)-min(dat$y)
                    peak.intensity.all <- max(peak.intensity.all, na.rm = TRUE)
                    peak.fit_1.contribution <- peak.intensity.fit_1/peak.intensity.all
                    peak.fit_2.contribution <- 1-peak.fit_1.contribution

                    # if(any(c(peak.fit_1.contribution*100, 
                    #          peak.fit_2.contribution*100) < 1)){
                    #      return(init.list[[4]] <- 2)
                    # } 
                    return(
                        list(
                            fit_1, fit_2, 
                            list("curve.one" = curve.one.contribution, 
                                "curve.two" = curve.two.contribution),
                            list("curve.one" = summary.fit.params["sigma1"],
                                "curve.two" = summary.fit.params["sigma2"]),
                            list("curve.one" = peak.fit_1.contribution,
                                 "curve.two" = peak.fit_2.contribution),                                
                            # list("curve.one" = summary.fit.params["C1"],
                            #      "curve.two" = summary.fit.params["C2"]),
                            init.list[length(init.list)][[1]]
                        )
                    )
                } else {
                    # fit curve 3
                    fit_3 <- draw.from.gaussian(
                        xs = x,
                        C = summary.fit.params["C3"],
                        SD = summary.fit.params["sigma3"]
                    )
                    peak.intensity.fit_3 <- max(fit_3, na.rm = TRUE)-min(dat$y, na.rm = TRUE)

                    # calculate area under each curve
                    integral.curve.three <- summary.fit.params["C3"]*summary.fit.params["sigma3"]*sqrt(2*pi)

                    # Calculate total area
                    integral.all.curves <- integral.curve.one+integral.curve.two+integral.curve.three

                    # percentage contribution of each curve
                    curve.one.contribution <- integral.curve.one/integral.all.curves
                    curve.two.contribution <- integral.curve.two/integral.all.curves
                    curve.three.contribution <- 1-curve.one.contribution-curve.two.contribution

                    # peak intensity contribution of each gaussian curve normalised
                    peak.intensity.all <- (peak.intensity.fit_1+peak.intensity.fit_2+peak.intensity.fit_3)
                    peak.fit_1.contribution <- peak.intensity.fit_1/peak.intensity.all
                    peak.fit_2.contribution <- peak.intensity.fit_2/peak.intensity.all
                    peak.fit_3.contribution <- 1-peak.fit_1.contribution-peak.fit_2.contribution

                    # print(paste(peak.fit_1.contribution*100, peak.fit_2.contribution*100, peak.fit_3.contribution*100))

                    if(any(c(peak.fit_1.contribution*100, 
                             peak.fit_2.contribution*100, 
                             peak.fit_3.contribution*100) < 5)){
                         return(init.list[[4]] <- 2)
                    }                    
                    return(
                        list(
                            fit_1, fit_2, fit_3, 
                            list("curve.one" = curve.one.contribution, 
                                "curve.two" = curve.two.contribution, 
                                "curve.three" = curve.three.contribution),
                            list("curve.one" = summary.fit.params["sigma1"],
                                "curve.two" = summary.fit.params["sigma2"],
                                "curve.three" = summary.fit.params["sigma3"]),
                            list("curve.one" = peak.fit_1.contribution,
                                 "curve.two" = peak.fit_2.contribution,
                                 "curve.three" = peak.fit_3.contribution),                                
                            # list("curve.one" = summary.fit.params["C1"],
                            #      "curve.two" = summary.fit.params["C2"],
                            #      "curve.three" = summary.fit.params["C3"]),                                
                            init.list[length(init.list)][[1]]
                        )
                    )
                }
            }
        },

        #' @description
        #' Make the plot of RMSD values with curves fitted to those values.
        #' @param dat Tibble of the RMSD values per k-mer.
        #' @param k Numeric vector of the k-mer.
        #' @param curve.vals Numeric vector of the values for the fitted curve to plot.
        #' @param nr.of.curves Numeric vector of the number of curves to fit.
        #' @param annot_plots if TRUE, will add 95% C.I. of each curve onto plot.
        #' @return list of plots.
        make_plot = function(dat, k, curve.vals, nr.of.curves, annot_plots){
            # extract standard deviation values to compute 95% confidence intervals
            st.devs <- unname(curve.vals[[length(curve.vals)-2]])
            CI.lst <- numeric(length = nr.of.curves*2)
            CI.lst[1] <- -1.96*st.devs[[1]]
            CI.lst[2] <- 1.96*st.devs[[1]]
            limits <- 500

            # if(any(CI.lst[2] > limits | CI.lst[1] < -limits)){
            #     return(2)
            # }

            if(nr.of.curves > 1){
                # if(any(1.96*st.devs[[2]] > limits | -1.96*st.devs[[2]] < -limits)){
                #     return(2)
                # }
                CI.lst[3] <- -1.96*st.devs[[2]]
                CI.lst[4] <- 1.96*st.devs[[2]]
            } 

            if(nr.of.curves == 3){
                # if(any(1.96*st.devs[[3]] > limits | -1.96*st.devs[[3]] < -limits)){
                #     return(2)
                # }
                CI.lst[5] <- -1.96*st.devs[[3]]
                CI.lst[6] <- 1.96*st.devs[[3]]
            }
            y.pos <- max(dat$y)

            # Plot the data with the model superimposed
            fit.plot <- dat %>%
                dplyr::mutate(kmer = 
                    stringr::str_replace_all(
                        string = stringr::str_to_title(kmer),
                        pattern = "_",
                        replacement = " "
                    )) %>% 
                ggplot(aes(x = x, y = y)) + 
                geom_line(alpha = 0.5) + 
                facet_wrap(~kmer) + 
                geom_line(
                    data = data.frame(
                        x = dat$x,
                        y = curve.vals[[1]]
                    ),
                    stat = "identity",
                    color = "#619B61",
                    alpha = 1,
                    linewidth = 2) + 
                theme_bw() + 
                theme(
                    axis.line = element_line(colour = "black"),
                    panel.background = element_blank(),
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    text = element_text(size = 25),
                    axis.title.y = element_blank()
                    # axis.text.y = element_blank(),
                    # axis.ticks.y = element_blank()
                ) + 
                labs(
                    x = "",
                    y = ""
                )

            if(annot_plots){
                fit.plot <- fit.plot + 
                    geom_vline(
                        xintercept = CI.lst[1], 
                        col = "#619B61", 
                        linetype = "dashed",
                        alpha = 0.7) +
                    geom_vline(
                        xintercept = CI.lst[2],
                        col = "#619B61", 
                        linetype = "dashed",
                        alpha = 0.7)
            } else {
                fit.plot <- fit.plot + 
                    coord_cartesian(xlim = c(-250, 250))
            }

            CI.annot <- function(nr.of.curves){
                if(nr.of.curves == 1){
                    return(paste0(
                        "95% CIs:\n",
                        "R: [", signif(CI.lst[1], 3), ",", 
                                    signif(CI.lst[2], 3), "]"
                    ))              
                } else if(nr.of.curves == 2){
                    return(paste0(
                        "95% CIs:\n",
                        "R: [", signif(CI.lst[1], 3), ",", 
                                    signif(CI.lst[2], 3), "]\n",
                        "B: [", signif(CI.lst[3], 3), ",", 
                                    signif(CI.lst[4], 3), "]"
                    ))              
                } else if(nr.of.curves == 3){
                    return(paste0(
                        "95% CIs:\n",
                        "R: [", signif(CI.lst[1], 3), ",", 
                                    signif(CI.lst[2], 3), "]\n",
                        "B: [", signif(CI.lst[3], 3), ",", 
                                    signif(CI.lst[4], 3), "]\n",
                        "G: [", signif(CI.lst[5], 3), ",", 
                                    signif(CI.lst[6], 3), "]"
                    ))
                }
            }

            if(nr.of.curves == 1){
                if(annot_plots){
                    fit.plot <- fit.plot + 
                        geom_text(
                            data = data.frame(
                                xpos = -Inf,
                                ypos = Inf,
                                annotateText = CI.annot(nr.of.curves = nr.of.curves),
                                hjustvar = -0.1, vjustvar = 1.1
                            ),
                            aes(
                                x = xpos,
                                y = ypos,
                                hjust = hjustvar,
                                vjust = vjustvar,
                                label = annotateText,
                                angle = 0
                            ),
                                size = 4
                            ) + 
                        theme(text = element_text(size = 25))
                } else {
                    fit.plot <- fit.plot + 
                            theme(
                                axis.line = element_line(colour = "black"),
                                strip.text.x = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                                panel.border = element_blank(),
                                text = element_text(size = 25),
                                axis.title.y = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank()
                            )
                }
                gaussian.linear.comb <- curve.vals[[1]]
            }

            if(nr.of.curves > 1){
                # Plot the data with the model superimposed
                fit.plot <- fit.plot + 
                    geom_line(
                    data = data.frame(
                        x = dat$x,
                        y = curve.vals[[2]]
                    ),
                    stat = "identity",
                    color = "#355c7d",
                    alpha = 1,
                    linewidth = 2)
                    # theme_bw() 
                    # theme(text = element_text(size = 25)) +
                    # labs(
                    #     x = "",
                    #     y = ""
                    # )

                if(annot_plots){
                    fit.plot <- fit.plot + 
                        geom_vline(
                            xintercept = CI.lst[3], 
                            col = "#355c7d", 
                            linetype = "dashed",
                            alpha = 0.7) +
                        geom_vline(
                            xintercept = CI.lst[4], 
                            col = "#355c7d", 
                            linetype = "dashed",
                            alpha = 0.7)
                }

                if(nr.of.curves == 2){
                    fit.plot <- fit.plot + 
                        geom_line(
                            data = data.frame(
                                x = dat$x, 
                                y = curve.vals[[1]]+curve.vals[[2]]-min(dat$y)
                            ),
                            stat = "identity",
                            color = "#f8b195",
                            alpha = 1,
                            linewidth = 2
                        )
                    
                    if(annot_plots){
                        fit.plot <- fit.plot + 
                            geom_text(
                                data = data.frame(
                                    xpos = -Inf,
                                    ypos = Inf,
                                annotateText = CI.annot(nr.of.curves = nr.of.curves),
                                hjustvar = -0.1, vjustvar = 1.1
                                ),
                                aes(
                                    x = xpos,
                                    y = ypos,
                                    hjust = hjustvar,
                                    vjust = vjustvar,
                                    label = annotateText,
                                    angle = 0
                                ),
                                    size = 4
                                )
                            # theme(text = element_text(size = 25))                        
                    } else {
                        fit.plot <- fit.plot + 
                            theme(
                                axis.line = element_line(colour = "black"),
                                strip.text.x = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                                panel.border = element_blank(),
                                text = element_text(size = 25),
                                axis.title.y = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank()                               
                            )
                    }
                    gaussian.linear.comb <- curve.vals[[1]]+curve.vals[[2]]-min(dat$y)
                } else {
                    fit.plot <- fit.plot + 
                        geom_line(
                            data = data.frame(
                                x = dat$x,
                                y = curve.vals[[3]]
                            ),
                            stat = "identity",
                            color = "#BB357E",
                            alpha = 1,
                            linewidth = 2
                        ) +
                        geom_line(
                            data = data.frame(
                                x = dat$x, 
                                y = curve.vals[[1]]+curve.vals[[2]]+curve.vals[[3]]-2*min(dat$y)
                            ),
                            stat = "identity",
                            color = "#f8b195",
                            alpha = 1,
                            linewidth = 2
                        ) 

                    if(annot_plots){
                        fit.plot <- fit.plot +
                            geom_vline(
                                xintercept = CI.lst[5], 
                                col = "#BB357E", 
                                linetype = "dashed",
                                alpha = 0.9) +
                            geom_vline(
                                xintercept = CI.lst[6], 
                                col = "#BB357E", 
                                linetype = "dashed",
                                alpha = 0.9) +
                            geom_text(
                                data = data.frame(
                                    xpos = -Inf,
                                    ypos = Inf,
                                    annotateText = CI.annot(nr.of.curves = nr.of.curves),
                                    hjustvar = -0.1, vjustvar = 1.1
                                ),
                                aes(
                                    x = xpos,
                                    y = ypos,
                                    hjust = hjustvar,
                                    vjust = vjustvar,
                                    label = annotateText,
                                    angle = 0
                                ),
                                size = 4
                            )
                            # theme(text = element_text(size = 25))
                    } else {
                        fit.plot <- fit.plot + 
                            theme(
                                axis.line = element_line(colour = "black"),
                                strip.text.x = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                                panel.border = element_blank(),
                                text = element_text(size = 25),
                                axis.title.y = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank()
                            )
                    }

                    gaussian.linear.comb <- curve.vals[[1]]+curve.vals[[2]]+curve.vals[[3]]-2*min(dat$y)
                }
            }
            return(list(fit.plot, CI.lst, gaussian.linear.comb))
        },

        #' @description
        #' Performs hierarchical clustering of RMSD tracts.
        #' @param kmer Numeric vector of the k-mer size.
        #' @return None.
        cluster_curves_into_ranges = function(kmer){
            t1 <- Sys.time()
            cur.msg <- paste0("Clustering RMSD tracts into ranges for kmer ", kmer)
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # get rmsd tracts
            all.files <- list.files(
                path = "../data",
                pattern = paste0("key_stats_kmer_", kmer, " *"), 
                recursive = TRUE,
                full.names = TRUE
            )
            all.curves <- c("curve.one", "curve.two", "curve.three")
            to.keep <- which(!grepl(pattern = "AverageRMSDFromClusters", x = all.files))
            all.files <- all.files[to.keep]

            results <- lapply(all.files, function(file){
                out <- fread(file)

                missing <- !(all.curves %in% colnames(out))
                if(any(missing)){
                    na.cols <- which(missing)
                    out <- cbind(out, rep(NA_real_, 4))
                    colnames(out)[length(colnames(out))] <- 
                        all.curves[na.cols[1]]

                    if(length(na.cols) > 1){
                        out <- cbind(out, rep(NA_real_, 4))
                        colnames(out)[length(colnames(out))] <- 
                            all.curves[na.cols[2]]
                    }
                }

                return(out)
            }) 
            results <- rbindlist(results)

            # Flatten all ranges; plot histogram/density plot
            df <- as_tibble(results) %>% 
                dplyr::filter(rowid == "ranges") %>% 
                dplyr::select(exp, rowid, curve.one, curve.two, curve.three) %>% 
                tidyr::gather(Curve, Value, -exp, -rowid) %>% 
                tidyr::drop_na(Value)

            #' @description 
            #' Log-transforms data for appropriate short/mid/long-range clustering.
            #' @param dat Tibble of the experiments with ranges.
            #' @param log_scale Boolean. If TRUE, log-transforms the data.
            #' @return None.
            cluster_ranges = function(dat, log_scale){
                dat$log_Value <- log2(dat$Value)

                # Clustering
                if(log_scale){
                    dat.dist <- dist(dat$log_Value)
                } else {
                    dat.dist <- dist(dat$Value)
                }

                clusts <- hclust(dat.dist, method = "ward.D2")
                cl_members <- cutree(clusts, k = 3)

                dat <- dat %>% 
                    dplyr::mutate(Cluster = as.numeric(cl_members))
                    
                hash_map <- dat %>% 
                    dplyr::group_by(Cluster) %>% 
                    dplyr::summarise(
                        lower.end = min(Value, na.rm = TRUE),
                        upper.end = max(Value, na.rm = TRUE)) %>% 
                    dplyr::mutate(Ranges = case_when(
                        upper.end == min(upper.end, na.rm = TRUE) ~ "short.range",
                        upper.end == median(upper.end, na.rm = TRUE) ~ "mid.range",
                        upper.end == max(upper.end, na.rm = TRUE) ~ "long.range")) %>% 
                    dplyr::pull(Ranges, Cluster)

                dat <- dat %>% 
                    dplyr::mutate(Cluster = case_when(
                        Cluster == names(hash_map)[1] ~ unname(hash_map)[1],
                        Cluster == names(hash_map)[2] ~ unname(hash_map)[2],
                        Cluster == names(hash_map)[3] ~ unname(hash_map)[3]
                    ))

                if(log_scale){
                    ranges.cluster.maps <- dat %>% 
                        dplyr::select(exp, Curve, Cluster) %>%
                        dplyr::arrange(exp, desc(Cluster))

                    # save curve/cluster map for gaussian curve re-mapping later
                    fwrite(
                        as.data.table(ranges.cluster.maps),
                        file = paste0("../data/ranges/kmer_", kmer, 
                                    "_Ranges_Mapped_to_Cluster.csv")
                    )

                    df <<- dat 
                }

                dat.plot <- dat %>% 
                    dplyr::mutate(
                        Cluster = dplyr::case_when(
                            Cluster == "short.range" ~ "Short range",
                            Cluster == "mid.range" ~ "Mid range",
                            Cluster == "long.range" ~ "Long range"
                        ),
                        hex = dplyr::case_when(
                            Cluster == "Short range" ~ "#619B61",
                            Cluster == "Mid range" ~ "#355B7D",
                            Cluster == "Long range" ~ "#BB357E"
                        )
                    ) %>%
                    ggplot(aes(x = Value)) +
                    geom_histogram(aes(
                        y = after_stat(density),
                        fill = hex, 
                        group = hex), 
                        colour = "black",
                        bins = 80
                    ) +
                    theme_bw() + 
                    theme_classic() + 
                    theme(text = element_text(size = 20)) + 
                    scale_fill_identity() +
                    labs(
                        title = paste(
                            "Sequence influences clustered into", 
                            length(unique(dat$Cluster)), 
                            "separate ranges"
                        ),
                        subtitle = ifelse(
                            log_scale, 
                            "Clustering based on log-scaled values", 
                            "Clustering based on true values"
                        ),
                        x = "Ranges (raw values)",
                        y = "Density in Cluster"
                    )
                    
                dir.create(
                    path = "../figures/ranges",
                    showWarnings = FALSE,
                    recursive = TRUE
                )
                ggsave(
                    plot = dat.plot,
                    filename = paste0("../figures/ranges/kmer_", kmer, 
                                        ifelse(log_scale, "_LOG_", "_"), 
                                        "flattened_ranges_histogram.pdf"),
                    width = 10, 
                    height = 7
                )

                # Hierarchical clustering plot
                pdf(paste0("../figures/ranges/kmer_", kmer, 
                            ifelse(log_scale, "_LOG_", "_"), 
                            "flattened_ranges_clustering.pdf"))
                plot(
                    clusts, 
                    hang = -1, 
                    cex = 0.2,
                    main = paste0("Hierarchical clustering on ", 
                                    "flattened ranges of sequence influences"),
                    xlab = "Ranges", 
                    sub = NA
                )

                rect.hclust(
                    clusts, 
                    k = max(cl_members, na.rm = TRUE), 
                    which = seq_len(max(cl_members, na.rm = TRUE)), 
                    border = seq_len(max(cl_members, na.rm = TRUE)), 
                    cluster = cl_members
                )
                save.plot <- dev.off()
            }

            # perform clustering and saves results
            for(scaling in c(TRUE, FALSE)){
                cluster_ranges(
                    dat = df,
                    log_scale = scaling
                )
            }

            # find cut-off values for each range
            cutoff.ranges <- df %>% 
                dplyr::group_by(Cluster) %>% 
                dplyr::summarise(
                    lower.end = min(Value, na.rm = TRUE),
                    upper.end = max(Value, na.rm = TRUE)
                ) %>% 
                dplyr::mutate(lower.end = ifelse(
                    Cluster == "short.range", upper.end, lower.end
                ))

            dir.create(
                path = "../data/ranges",
                showWarnings = FALSE,
                recursive = TRUE
            )
            fwrite(
                x = cutoff.ranges, 
                file = paste0("../data/ranges/kmer_", kmer, 
                                "_Ranges_cutoffs_from_clustering.csv")
            )

            # get ranges for each exp
            ranges.df <- df %>% 
                dplyr::group_by(Cluster) %>% 
                dplyr::select(-c(Curve, log_Value))

            ranges.df <- ranges.df %>% 
                dplyr::mutate(exp = factor(exp, levels = c(
                    base::unique(ranges.df$exp) %>% 
                        stringr::str_sort(., numeric = TRUE)
                ))) %>% 
                dplyr::arrange(exp, Value) %>% 
                dplyr::group_by(exp, Cluster) %>% 
                dplyr::mutate(
                    Cluster = ifelse(Value == min(Value, na.rm = TRUE), 
                        paste(Cluster, "-1", sep = ""), 
                        paste(Cluster, "-2", sep = ""))) %>% 
                dplyr::ungroup() %>% 
                tidyr::unite(
                    col = "ID",
                    c(exp, rowid, Value),
                    sep = "_",
                    remove = FALSE
                ) 
            
            temp <- ranges.df %>% 
                dplyr::select(ID, Cluster)

            join_tables <- function(metric){                
                other.metric <- results %>% 
                    dplyr::filter(rowid == "ranges" | rowid == {{metric}}) %>% 
                    dplyr::select(exp, rowid, curve.one, curve.two, curve.three) %>% 
                    tidyr::gather(Curve, Value, -exp, -rowid) %>% 
                    tidyr::drop_na(Value) %>% 
                    dplyr::group_by(exp) %>%
                    dplyr::mutate(exp = factor(exp, 
                        levels = levels(ranges.df$exp)
                    )) %>% 
                    tidyr::unite(
                        col = "ID",
                        c(exp, rowid, Value),
                        sep = "_",
                        remove = FALSE
                    ) 
                    
                other.metric <- dplyr::left_join(
                        x = other.metric,
                        y = temp,
                        by = "ID"
                    ) %>% 
                    dplyr::group_by(exp) %>% 
                    dplyr::mutate(Cluster = 
                        ifelse(is.na(Cluster) & metric == "peak.intensity", 
                            lag(Cluster), 
                        ifelse(is.na(Cluster) & metric == "contribution", 
                            lead(Cluster), Cluster
                        ))
                    ) %>%
                    dplyr::select(-ID) %>% 
                    dplyr::mutate(
                        Cluster = factor(Cluster, levels = c(
                            "short.range-1", "short.range-2",
                            "mid.range-1", "mid.range-2",
                            "long.range-1", "long.range-2"
                        ))
                    ) %>%
                    dplyr::arrange(exp, Cluster) 

                return(other.metric)
            }
            df.pi <- join_tables(metric = "peak.intensity")
            df.c <- join_tables(metric = "contribution")

            # reorder elements
            ranges.df <- df.pi %>% 
                dplyr::bind_rows(df.c) %>% 
                dplyr::group_by(exp, Cluster) %>% 
                dplyr::arrange(exp, desc(rowid)) %>%
                dplyr::distinct()

            ranges.final.df <- ranges.df %>%  
                tidyr::spread(Cluster, Value) %>% 
                dplyr::select(
                    stringr::str_sort(names(.), numeric = TRUE)) %>% 
                dplyr::ungroup()            

            fwrite(
                x = ranges.final.df, 
                file = paste0("../data/ranges/kmer_", kmer, 
                            "_Ranges_cutoffs_from_clustering", 
                            "_all-exp.csv")
            )

            files.to.import <- list.files(
                path = "../data/ranges", 
                pattern = "_Ranges_cutoffs_from_clustering.csv",
                full.names = TRUE
            )

            kmer.val <- as.integer(stringr::str_extract(
                string = basename(files.to.import), 
                pattern = "[:digit:]"
            ))

            if(kmer == 8){
                out <- lapply(1:length(files.to.import), function(i){
                    out <- fread(files.to.import[i])
                    out <- cbind(out, "kmer" = kmer.val[i])
                })
                out <- rbindlist(out)
                out[, lower.end := NULL]

                out <- as_tibble(out) %>% 
                    dplyr::mutate(kmer = paste0("kmer_", kmer)) %>% 
                    tidyr::spread(key = "kmer", value = "upper.end")

                fwrite(
                    x = out, 
                    file = "../data/ranges/Ranges_cutoffs.csv"
                )
            }

            # count number of curves to fit per exp
            curve.counts <- df %>% 
                dplyr::group_by(exp) %>%
                dplyr::summarise(Count = dplyr::n())

            fwrite(
                x = curve.counts, 
                file = paste0("../data/ranges/kmer_", kmer, 
                                "_curve_counts.csv")
            )

            if(kmer == 8){
                # Reassign org_file with curve counts per kmer
                files.to.import <- list.files(
                    path = "../data/ranges", 
                    pattern = "curve_counts",
                    full.names = TRUE
                )
                out <- lapply(1:length(files.to.import), function(i){
                    out <- fread(files.to.import[i])
                    out <- cbind(out, "kmer" = kmer.val[i])
                })
                out <- rbindlist(out)

                out <- as_tibble(out) %>% 
                    dplyr::mutate(kmer = paste0("kmer_", kmer)) %>% 
                    tidyr::spread(key = "kmer", value = "Count")

                # improt org_file.csv
                org.file <- fread("../../data/org_file.csv")
                kmer.cols <- which(grepl("^kmer_", colnames(org.file)))
                if(any(kmer.cols)){
                    org.file[, which(grepl("^kmer_", colnames(org.file))) := NULL]
                    fwrite(x = org.file, file = "../../data/org_file.csv")
                }
                system("cp ../../data/org_file.csv ../../data/org_file_backup.csv")

                org.file <- as_tibble(org.file) %>% 
                    dplyr::mutate(
                        exp = paste0(`Fragmentation type`, "/", `Experiment folder`), 
                        .before = 1
                    ) %>% 
                    dplyr::left_join(., out, by = "exp")

                org.file <- org.file %>% 
                    dplyr::select(-exp)

                fwrite(
                    x = org.file, 
                    file = "../../data/org_file.csv"
                )
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Performs hierarchical clusterings of RMSD tracts.
        #' @param kmer Numeric vector of the k-mer size.
        #' @return None.
        get_rmsd_tracts = function(kmer){
            all.files <- list.files(
                path = "../data",
                pattern = paste0("key_stats_kmer_", kmer, ".*"), 
                recursive = TRUE,
                full.names = TRUE
            )

            if(private$auto_fit){
                curve.labels <- c("long.range", "mid.range", "short.range")
            } else {
                curve.labels <- c(
                    "short.range-1", "short.range-2",
                    "mid.range-1", "mid.range-2",
                    "long.range-1", "long.range-2"
                )
            }

            results <- lapply(all.files, function(file){
                out <- fread(file, check.names = TRUE)

                if(private$auto_fit){
                    if(!("curve.two" %in% colnames(out))){
                        out <- cbind(out, rep(NA_real_, 4))
                        colnames(out)[length(colnames(out))] <- "curve.two"
                    }

                    if(!("curve.three" %in% colnames(out))){
                        out <- cbind(out, rep(NA_real_, 4))
                        colnames(out)[length(colnames(out))] <- "curve.three"
                    }
                } else {
                    not.duplicate.col.names <- str_extract(
                        string = colnames(out)[5:length(out)], 
                        pattern = "[:digit:]"
                    )

                    colnames(out)[5:length(out)][which(is.na(not.duplicate.col.names))] <- 
                        paste0(colnames(out)[5:length(out)][which(is.na(not.duplicate.col.names))], "-1")

                    colnames(out)[5:length(out)][which(!is.na(not.duplicate.col.names))] <- 
                    stringr::str_replace(
                        string = colnames(out)[5:length(out)][which(!is.na(not.duplicate.col.names))], 
                        pattern = "\\.[0-9]",
                        replacement = "-2"
                    )

                    missing <- which(!curve.labels %in% colnames(out))
                    if(any(missing)){
                        for(m in missing){
                            out <- cbind(out, rep(NA_real_, 4))
                            colnames(out)[length(colnames(out))] <- curve.labels[m]
                        }
                    }
                }
                return(out)
            })
            results <- do.call(rbind, results)
            results <- as_tibble(results)

            if(private$auto_fit){
                results <- results %>% 
                    dplyr::relocate(
                        curve.one, 
                        curve.two, 
                        curve.three, 
                        .after = 4) %>% 
                    dplyr::arrange(category)

                df <- results %>% 
                    dplyr::filter(rowid == "ranges") %>% 
                    dplyr::select(-c(
                        category, 
                        category_general, 
                        rowid)) %>% 
                    dplyr::relocate(
                        curve.one, 
                        curve.two, 
                        curve.three, 
                        .after = 1
                    )

                colmsn <- colMeans(df[, -1], na.rm = TRUE)
                df <- as_tibble(cbind(df[, "exp"], df[, -1][order(colmsn, decreasing = FALSE)]))
            } else {
                rearrange.cols <- 4+match(colnames(results[, 5:length(results)]), curve.labels)

                results <- results %>% 
                    dplyr::select(order(c(
                        1, 2, 3, 4, 
                        rearrange.cols))) %>% 
                    dplyr::select(where(~sum(!is.na(.)) > 0))

                df <- results %>% 
                    dplyr::filter(rowid == "ranges") %>% 
                    dplyr::select(-c(
                        category, 
                        category.sub, 
                        rowid)
                    )
            }

            dir.create(
                path = paste0("../data/seq_influence"),
                showWarnings = FALSE
            )
            fwrite(
                results,
                paste0("../data/seq_influence/", 
                        ifelse(private$auto_fit, "key", "new"), 
                        "_stats_kmer_", kmer, ".csv")
            )

            # hierarchical clustering
            df.hc <- apply(df[,-1], 2, as.numeric)
            rownames(df.hc) <- df$exp
            df.dist <- dist(df.hc) %>% suppressWarnings()
            df.dendro <- as.dendrogram(hclust(df.dist, method = "ward.D2"))
            transpose.df <- t(df.hc)

            dir.create(
                path = "../figures/seq_influence", 
                showWarnings = FALSE
            )
            pdf(
                paste0("../figures/seq_influence/", kmer, "-mer-clustering_", 
                       ifelse(private$auto_fit, "auto_fit", "optimised"),
                       ".pdf"), 
                height = 10, 
                width = 8
            )

            cols <- 200
            heatmap.breaks <- seq(min(df.hc, na.rm = TRUE), 
                                  max(df.hc, na.rm = TRUE), 
                                  length.out = cols+1)
            heatmap.2(
                df.hc,
                offsetRow = 0,
                offsetCol = 0,
                Rowv = df.dendro,
                Colv = FALSE,
                dendrogram = "row",
                revC = FALSE,
                trace = "none",
                density.info = "histogram",
                col = hcl.colors(cols, "RdYlGn"), 
                breaks = heatmap.breaks,
                na.color = "grey",
                notecol = "black",
                cexCol = 0.4,
                cexRow = 0.4,
                labRow = rownames(df.hc),
                labCol = colnames(df.hc),
                margins = c(4,20),
                key.xlab = "RMSD ranges"
            )
            plot.save <- dev.off()

            if(!private$auto_fit){
                # Plot individual distributions for each range
                p <- results %>% 
                    dplyr::filter(rowid == "ranges") %>% 
                    dplyr::select(
                        exp, 
                        category, 
                        contains("range")) %>% 
                    tidyr::gather(
                        -c(exp,category), 
                        key = "Ranges",
                        value = "Value") %>% 
                    dplyr::mutate(
                        Ranges = factor(Ranges, 
                        levels = c("short.range", 
                                   "mid.range", 
                                   "long.range"))) %>% 
                    ggplot(aes(x = Value)) + 
                    geom_histogram(aes(
                        fill = category),
                        bins = 50,
                        colour = "black"
                    ) +
                    facet_wrap(~Ranges, scales = "free_x") +
                    labs(
                        title = "Sequence influences clustered into 3 separate ranges",
                        subtitle = paste0("Kmer: ", kmer),
                        x = "Ranges",
                        y = "Count"
                    )

                suppressWarnings(ggsave(
                    plot = p, 
                    filename = paste0("../figures/seq_influence/", 
                                      kmer, "-mer-ranges.pdf"),
                    height = 7, 
                    width = 13
                ))
            }

            if(!private$auto_fit){
                # categorical hierarchical clustering
                df <- df %>% 
                    tidyr::gather(-exp, key = "Curve", value = "Value") %>% 
                    dplyr::mutate(Value = case_when(
                        (str_detect(Curve, "short") & !is.na(Value)) ~ 1,
                        (str_detect(Curve, "mid") & !is.na(Value)) ~ 2,
                        (str_detect(Curve, "long") & !is.na(Value)) ~ 3)) %>% 
                    tidyr::spread(Curve, Value)

                rearrange.cols <- match(colnames(df[-1]), curve.labels)

                df <- df %>% 
                    dplyr::select(order(c(1, rearrange.cols))) %>% 
                    dplyr::select(where(~sum(!is.na(.)) > 0))
                
                df.hc <- df[-1] %>% 
                    dplyr::mutate(dplyr::across(
                        where(is.numeric), 
                        function(x) ifelse(is.na(x), 0, x))) %>% 
                    dplyr::mutate(dplyr::across(
                        where(is.numeric), 
                        as.factor)
                    )

                df.dist <- daisy(df.hc, metric = "gower") %>% suppressWarnings()
                df.dendro <- as.dendrogram(hclust(df.dist, method = "ward.D2"))
                df.hc <- apply(df.hc, 2, as.numeric)
                rownames(df.hc) <- df$exp

                pdf(
                    paste0("../figures/seq_influence/", kmer, 
                           "-mer-clustering_optimised-categorical.pdf"), 
                    height = 10, 
                    width = 8
                )

                breaks <- -1:3
                col <- c("grey","red","orange", "darkgreen")

                heatmap.2(
                    df.hc,
                    offsetRow = 0,
                    offsetCol = 0,
                    Rowv = df.dendro,
                    Colv = FALSE,
                    dendrogram = "row",
                    revC = FALSE,
                    trace = "none",
                    density.info = "histogram",
                    col = col, 
                    breaks = breaks,
                    na.color = "grey",
                    notecol = "black",
                    cexCol = 0.4,
                    cexRow = 0.4,
                    labRow = rownames(df.hc),
                    labCol = colnames(df.hc),
                    margins = c(4,20),
                    key.xlab = "Category of curve"
                )
                plot.save <- dev.off()
            }
        }
    )
)