FitCurves <- R6::R6Class(
    classname = "FitCurves",
    public = list(
        #' @field chr Numeric vector of chromosome number.
        chr = NULL,

        #' @field which_exp_ind Numeric vector of the index of the experiment 
        #' to process from org_file. If NULL, all exp will be processed.
        which_exp_ind = NULL,

        initialize = function(chr, which_exp_ind){
            if(!missing(chr)) self$chr <- chr
            if(!missing(which_exp_ind)) self$which_exp_ind <- which_exp_ind
            
            # get full org_file.csv
            private$get_org_file()

            if(is.null(self$which_exp_ind)){
                self$which_exp_ind <- which(
                    (private$org_file$`DSB Map` == "TRUE") & 
                    (private$org_file$`RMSD?` == "TRUE")
                )
            }
            private$len_of_loop <- self$which_exp_ind
        },

        #' @description
        #' Generate plots from RMSD calculations for each experiment.
        #' @return None.
        generate_rmsd_plots = function(){
            start.time <- Sys.time()
            cur.msg <- "Generating RMSD plots and fitting Gaussian curves"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")

            for(i in private$len_of_loop){
                private$bp_exp <- paste0(
                    private$org_file[i, `Fragmentation type`], "/",
                    private$org_file[i, `Experiment folder`]
                )
                private$category <- private$org_file[i, Category_Main]
                private$category_sub <- private$org_file[i, Category_sub]

                cur.msg <- paste0("Processing ", i, "/", 
                                  nrow(private$org_file),
                                  ". ", private$bp_exp)
                cat(cur.msg, "\n", sep = "")

                private$fits <- c("kmer_4" = 3, "kmer_6" = 3, "kmer_8" = 3)
                private$plot_rmsd()
            }

            if(is.null(self$which_exp_ind)){
                for(kmer in c(4,6,8)){                    
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

        #' @field category_sub Character vector of the breakage sub-category.
        category_sub = NULL,

        #' @field auto_fit Boolean. If TRUE, will not optimise 
        #' number of curves to fit.
        auto_fit = TRUE,

        #' @field len_of_loop Numeric vector of how many experiments to process.
        len_of_loop = NULL,

        #' @field fits Named vector of curves to fit per k-mer RMSD values.
        fits = NULL,

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
        #' Generate RMSD plots for all k-mers available
        #' @return None.
        plot_rmsd = function(){
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
                ggplot(aes(x = x, y = y)) + 
                    geom_line(linewidth = 0.8) + 
                    facet_wrap(~kmer, ncol = 4, scales = "free_y") + 
                    theme_bw() + 
                    coord_cartesian(ylim = c(0, NA)) + 
                    labs(
                        x = "Position away from breakpoint",
                        y = "RMSD"
                    )

            dir.create(
                path = paste0("../figures/", private$bp_exp),
                showWarnings = FALSE,
                recursive = TRUE
            )
            ggsave(
                filename = paste0("../figures/", private$bp_exp, 
                                    "/chr", self$chr, 
                                    "_sidebyside_RMSD.pdf"),
                plot = plots, 
                height = 6,
                width = 21
            )

            # GMM plots
            data.sets <- data.sets %>%
                dplyr::mutate(dplyr::across(where(is.factor), as.character)) %>% 
                dplyr::filter(
                    kmer == "kmer_2" |
                    kmer == "kmer_4" |
                    kmer == "kmer_6" |
                    kmer == "kmer_8"
                )

            p <- vector(mode = "list", length = 3)
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

                    output <- private$make_plot(
                        dat = dat, 
                        k = names(private$fits[k]), 
                        curve.vals = curvefits, 
                        nr.of.curves = private$fits[k]
                    )
                    p[[k]] <- output[[1]]

                    # percent contribution of each gaussian curve towards breakage
                    df <- as_tibble(curvefits[length(curvefits)-2][[1]])
                    df[1,] <- df[1,]*100

                    # range of influence based on 95 percent confidence intervals
                    output.CIlst <- output[[2]]
                    curves <- output.CIlst[seq(2, length(output.CIlst), 2)]-
                              output.CIlst[seq(1, length(output.CIlst), 2)]

                    # sd of curves
                    SD <- as.numeric(curvefits[[length(curvefits)-1]])
                    SD <- c(SD, rep(NA_real_, 3-length(SD)))

                    # peak intensity of each gaussian curve
                    peak.intensity <- sapply(1:private$fits[k], function(x){
                        max(curvefits[[x]])
                    })
                    peak.intensity <- c(
                        peak.intensity, 
                        rep(NA_real_, 3-length(peak.intensity))
                    )

                    # combine results
                    df <- rbind(df, curves, SD, peak.intensity)
                    df <- df %>% 
                        dplyr::mutate(
                            exp = private$bp_exp,
                            category = private$category,
                            category_sub = private$category_sub,
                            rowid = c("contribution", "ranges", "SD", "peak.intensty"), 
                            .before = 1,
                        )

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

                    dir.create(
                        path = paste0("../data/", private$bp_exp),
                        showWarnings = FALSE,
                        recursive = TRUE
                    )
                    fwrite(
                        df,
                        file = paste0("../data/", private$bp_exp, "/key",
                                      "_stats_", names(private$fits[k]), ".csv"),
                        row.names = FALSE
                    )
                }
                # if(k == 1) print(p) # temp. avoid an unsolved ggplot bug
            }

            ggsave(
                filename = paste0("../figures/", private$bp_exp, 
                                  "/chr", self$chr, "_RMSD_all_GMMFit.pdf"),
                plot = cowplot::plot_grid(plotlist = p, axis = "b", ncol = 1), 
                height = 15, 
                width = 7
            )
        },

        #' @description
        #' Fit gaussian curves to underlying RMSD values.
        #' @param dat Tibble of the RMSD values per k-mer.
        #' @param ind Numeric vector of the index.
        #' @param C.value Numeric vector of the starting gaussian coefficient.
        #' @param sigma3 Numeric vector of the starting sigma coefficient.
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
                sd1 = seq(2, 16, 0.5),
                sd2 = seq(8, 50, 1.5),
                sd3 = seq(20, 77, 2),
                rss = rep(NA, length(seq(2, 16, 0.5)))
            )
            
            nls.res <- lapply(1:nrow(grid_search), function(index){
                if(nr.of.curves == 2){
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
                        upper=rep(c(NULL, 248), 3),
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
                        upper=rep(c(NULL, 248), 3),
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

            # # check for zero coefficients
            # any.zero.coefs <- any(summary.fit.params[grepl(
            #     pattern = "C", x = names(summary.fit.params)
            # )] == 0)
            # if(any.zero.coefs) return(init.list[[4]] <- 2)
            
            # fit each Gaussian curve separately to data
            draw.from.gaussian <- function(xs, C, SD, miny){
                return(C*exp(-(xs-0)^2/(2*SD^2))+min(miny))
            }

            # integral functions for curve 1
            curve.one <- function(x, C, SD){
                return(
                    (C*exp(-(x-0)^2/(2*SD^2)))
                )
            }

            # integral functions for curve 2
            curve.two <- function(x, C, SD){
                return(
                    (C*exp(-(x-0)^2/(2*SD^2)))
                )
            }

            # integral functions for curve 3
            curve.three <- function(x, C, SD){
                return(
                    (C*exp(-(x-0)^2/(2*SD^2)))
                )
            }

            # integral functions for all curves
            all.curves <- function(x, nr.of.curves, C1, C2, 
                                   C3=NULL, SD1, SD2, SD3=NULL){
                if(nr.of.curves == 2){
                    return(
                        (C1*exp(-(x-0)^2/(2*SD1^2))+
                        C2*exp(-(x-0)^2/(2*SD2^2))
                        )
                    )
                } else if(nr.of.curves == 3){
                    return(
                        (C1*exp(-(x-0)^2/(2*SD1^2))+
                        C2*exp(-(x-0)^2/(2*SD2^2))+
                        C3*exp(-(x-0)^2/(2*SD3^2))
                        )
                    )
                }
            }
            
            # fit curve 1
            fit_1 <- draw.from.gaussian(
                xs = x,
                C = summary.fit.params["C1"],
                SD = summary.fit.params["sigma1"],
                miny = min(y, na.rm = TRUE)
            )

            if(nr.of.curves == 1){
                return(
                    list(
                        fit_1, 
                        list("curve.one" = 1),
                        list("curve.one" = summary.fit.params["sigma1"]),
                        init.list[length(init.list)][[1]]
                    )
                )
            } else if(nr.of.curves > 1){
                # fit curve 2
                fit_2 <- draw.from.gaussian(
                    xs = x,
                    C = summary.fit.params["C2"],
                    SD = summary.fit.params["sigma2"],
                    miny = min(y, na.rm = TRUE)
                )

                # integrate curve 1
                integral.curve.one <- tryCatch({
                    integrate(
                        curve.one, 
                        lower = x[1], 
                        upper = x[length(x)],
                        C = summary.fit.params["C1"],
                        SD = summary.fit.params["sigma1"])$value
                    }, error = function(e){
                    integrate(
                        curve.one, 
                        lower = x[1], 
                        upper = x[length(x)],
                        C = summary.fit.params["C1"],
                        SD = summary.fit.params["sigma1"],
                        rel.tol = 1e-15)$value
                })

                # integrate curve 2
                integral.curve.two <- tryCatch({
                    integrate(
                        curve.two, 
                        lower = x[1], 
                        upper = x[length(x)],
                        C = summary.fit.params["C2"],
                        SD = summary.fit.params["sigma2"])$value
                    }, error = function(e){
                    integrate(
                        curve.two, 
                        lower = x[1], 
                        upper = x[length(x)],
                        C = summary.fit.params["C2"],
                        SD = summary.fit.params["sigma2"],
                        rel.tol = 1e-15)$value
                })

                if(nr.of.curves == 2){
                    # integrate linear combination of gaussian curves
                    integral.all.curves <- tryCatch({
                        integrate(
                            all.curves, 
                            lower = x[1], 
                            upper = x[length(x)],
                            nr.of.curves = 2,
                            C1 = summary.fit.params["C1"],
                            C2 = summary.fit.params["C2"],
                            SD1 = summary.fit.params["sigma1"],
                            SD2 = summary.fit.params["sigma2"])$value
                    },error = function(e){
                        integrate(
                            all.curves, 
                            lower = x[1], 
                            upper = x[length(x)],
                            nr.of.curves = 2,
                            C1 = summary.fit.params["C1"],
                            C2 = summary.fit.params["C2"],
                            SD1 = summary.fit.params["sigma1"],
                            SD2 = summary.fit.params["sigma2"],
                            rel.tol = 1e-15)$value
                    })

                    # percentage contribution of each curve
                    curve.one.contribution <- integral.curve.one/integral.all.curves
                    curve.two.contribution <- 1-curve.one.contribution
                    return(
                        list(
                            fit_1, fit_2, 
                            list("curve.one" = curve.one.contribution, 
                                "curve.two" = curve.two.contribution),
                            list("curve.one" = summary.fit.params["sigma1"],
                                "curve.two" = summary.fit.params["sigma2"]),
                            init.list[length(init.list)][[1]]
                        )
                    )
                } else {
                    # fit curve 3
                    fit_3 <- draw.from.gaussian(
                        xs = x,
                        C = summary.fit.params["C3"],
                        SD = summary.fit.params["sigma3"],
                        miny = min(y, na.rm = TRUE)
                    )

                    # integrate curve 3
                    integral.curve.three <- tryCatch({
                        integrate(
                            curve.three, 
                            lower = x[1], 
                            upper = x[length(x)],
                            C = summary.fit.params["C3"],
                            SD = summary.fit.params["sigma3"])$value
                    }, error = function(e){
                        integrate(
                            curve.three, 
                            lower = x[1], 
                            upper = x[length(x)],
                            C = summary.fit.params["C3"],
                            SD = summary.fit.params["sigma3"],
                            rel.tol = 1e-15)$value
                    })

                    # integrate linear combination of gaussian curves
                    integral.all.curves <- tryCatch({
                        integrate(
                            all.curves, 
                            lower = x[1], 
                            upper = x[length(x)],
                            nr.of.curves = 3,
                            C1 = summary.fit.params["C1"],
                            C2 = summary.fit.params["C2"],
                            C3 = summary.fit.params["C3"],
                            SD1 = summary.fit.params["sigma1"],
                            SD2 = summary.fit.params["sigma2"],
                            SD3 = summary.fit.params["sigma3"])$value
                    }, error = function(e){
                        integrate(
                            all.curves, 
                            lower = x[1], 
                            upper = x[length(x)],
                            nr.of.curves = 3,
                            C1 = summary.fit.params["C1"],
                            C2 = summary.fit.params["C2"],
                            C3 = summary.fit.params["C3"],
                            SD1 = summary.fit.params["sigma1"],
                            SD2 = summary.fit.params["sigma2"],
                            SD3 = summary.fit.params["sigma3"],
                            rel.tol = 1e-15
                            )$value
                    })

                    # percentage contribution of each curve
                    curve.one.contribution <- integral.curve.one/integral.all.curves
                    curve.two.contribution <- integral.curve.two/integral.all.curves
                    curve.three.contribution <- 1-curve.one.contribution-curve.two.contribution

                    if(any(c(curve.one.contribution, 
                             curve.two.contribution, 
                             curve.three.contribution) < 0.05)){
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
        #' @return list of plots.
        make_plot = function(dat, k, curve.vals, nr.of.curves){
            # extract standard deviation values to compute 95% confidence intervals
            st.devs <- unname(curve.vals[[length(curve.vals)-1]])
            CI.lst <- numeric(length = nr.of.curves*2)
            CI.lst[1] <- -1.96*st.devs[[1]]
            CI.lst[2] <- 1.96*st.devs[[1]]
            limits <- 500

            if(any(CI.lst[2] > limits | CI.lst[1] < -limits)){
                return(2)
            }

            if(nr.of.curves > 1){
                if(any(1.96*st.devs[[2]] > limits | -1.96*st.devs[[2]] < -limits)){
                    return(2)
                }
                CI.lst[3] <- -1.96*st.devs[[2]]
                CI.lst[4] <- 1.96*st.devs[[2]]
            } 

            if(nr.of.curves == 3){
                if(any(1.96*st.devs[[3]] > limits | -1.96*st.devs[[3]] < -limits)){
                    return(2)
                }
                CI.lst[5] <- -1.96*st.devs[[3]]
                CI.lst[6] <- 1.96*st.devs[[3]]
            }
            y.pos <- max(dat$y)

            # Plot the data with the model superimposed
            fit.plot <- dat %>%
                ggplot(aes(x = x, y = y)) + 
                geom_line(alpha = 0.8) + 
                facet_wrap(~kmer) + 
                geom_line(
                    data = data.frame(
                        x = dat$x,
                        y = curve.vals[[1]]
                    ),
                    stat = "identity",
                    color = "red",
                    alpha = 0.8,
                    linewidth = 1.2) + 
                geom_vline(
                    xintercept = CI.lst[1], 
                    col = "red", 
                    linetype = "dashed") + 
                geom_vline(
                    xintercept = CI.lst[2],
                    col = "red", 
                    linetype = "dashed") +
                theme_bw() + 
                theme(
                    strip.text.x = element_text(size = 10),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 10)
                ) + 
                labs(
                    x = "",
                    y = ""
                )

            CI.annot <- function(nr.of.curves){
                if(nr.of.curves == 2){
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

            # if(nr.of.curves == 1){
            #     fit.plot <- fit.plot + 
            #     geom_text(
            #         data = data.frame(
            #             xpos = -Inf,
            #             ypos = Inf,
            #             annotateText = paste0(
            #                 "Red CI:   ", formatC(signif(CI.lst[1], 3), 3), 
            #                 " / ", formatC(signif(CI.lst[2], 3), 3)
            #             ),
            #             hjustvar = -0.1, vjustvar = 1.1
            #         ),
            #         aes(
            #             x = xpos,
            #             y = ypos,
            #             hjust = hjustvar,
            #             vjust = vjustvar,
            #             label = annotateText,
            #             angle = 0
            #         ),
            #         size = 4
            #     )
            # }

            if(nr.of.curves > 1){
                # Plot the data with the model superimposed
                fit.plot <- fit.plot + 
                    geom_line(
                    data = data.frame(
                        x = dat$x,
                        y = curve.vals[[2]]
                    ),
                    stat = "identity",
                    color = "blue",
                    alpha = 0.8,
                    linewidth = 1.2) +
                    geom_vline(
                        xintercept = CI.lst[3], 
                        col = "blue", 
                        linetype = "dashed") + 
                    geom_vline(
                        xintercept = CI.lst[4], 
                        col = "blue", 
                        linetype = "dashed") + 
                    theme_bw() + 
                    theme(
                        strip.text.x = element_text(size = 10),
                        axis.text.x = element_text(size = 10),
                        axis.text.y = element_text(size = 10)
                    ) + 
                    labs(
                        x = "",
                        y = ""
                    )
                if(nr.of.curves == 2){
                    fit.plot <- fit.plot + 
                        geom_line(
                            data = data.frame(
                                x = dat$x, 
                                y = curve.vals[[1]]+curve.vals[[2]]-min(dat$y)
                            ),
                                stat = "identity",
                                color = "orange",
                                alpha = 0.8,
                                linewidth = 1) + 
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
                } else {
                    fit.plot <- fit.plot + 
                        geom_vline(
                            xintercept = CI.lst[5], 
                            col = "darkgreen", 
                            linetype = "dashed") + 
                        geom_vline(
                            xintercept = CI.lst[6], 
                            col = "darkgreen", 
                            linetype = "dashed"
                        ) +
                        geom_line(
                            data = data.frame(
                                x = dat$x,
                                y = curve.vals[[3]]
                            ),
                            stat = "identity",
                            color = "darkgreen",
                            alpha = 0.8,
                            linewidth = 1.2) + 
                        geom_line(
                            data = data.frame(
                                x = dat$x, 
                                y = curve.vals[[1]]+curve.vals[[2]]+curve.vals[[3]]-2*min(dat$y)
                            ),
                            stat = "identity",
                            color = "orange",
                            alpha = 0.8,
                            linewidth = 1) + 
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
                }
            }
            return(list(fit.plot, CI.lst))
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

            results <- lapply(all.files, function(file){
                out <- fread(file)
                out <- as_tibble(out)

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
            results <- do.call(rbind, results)

            # Flatten all ranges; plot histogram/density plot
            df <- results %>% 
                dplyr::filter(rowid == "ranges") %>% 
                dplyr::select(exp, curve.one, curve.two, curve.three) %>% 
                tidyr::gather(-exp, key = "Curve", value = "Value") %>% 
                dplyr::filter(!is.na(Value))

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

                clusts <- hclust(dat.dist, method = "complete")
                cl_members <- cutree(clusts, k = 3)

                dat <- dat %>% 
                    dplyr::mutate(Cluster = as.numeric(cl_members))
                    
                hash_map <- dat %>% 
                    dplyr::group_by(Cluster) %>% 
                    dplyr::summarise(
                        lower.end = min(Value),
                        upper.end = max(Value)
                    ) %>% 
                    dplyr::mutate(Ranges = case_when(
                        upper.end == min(upper.end) ~ "short.range",
                        upper.end == median(upper.end) ~ "mid.range",
                        upper.end == max(upper.end) ~ "long.range",
                    )) %>% 
                    dplyr::pull(Ranges, Cluster)

                dat <- dat %>% 
                    dplyr::mutate(Cluster = case_when(
                        Cluster == names(hash_map)[1] ~ unname(hash_map)[1],
                        Cluster == names(hash_map)[2] ~ unname(hash_map)[2],
                        Cluster == names(hash_map)[3] ~ unname(hash_map)[3]
                    ))
                
                if(log_scale) df <<- dat 

                dat.plot <- dat %>% 
                    ggplot(aes(x = Value)) +
                    geom_histogram(aes(
                        y = after_stat(density),
                        fill = Cluster, 
                        group = Cluster), 
                        colour = "black",
                        bins = 130
                    ) +
                    # scale_fill_discrete(labels = label) +   
                    labs(
                        title = paste("Sequence influences clustered into", 
                                        length(unique(dat$Cluster)), "separate ranges"),
                        subtitle = ifelse(log_scale, 
                                            "Clustering based on log-scaled values", 
                                            "Clustering based on true values"),
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

                if(!log_scale){
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
                        k = max(cl_members), 
                        which = seq_len(max(cl_members)), 
                        border = seq_len(max(cl_members)), 
                        cluster = cl_members
                    )
                    save.plot <- dev.off()
                }
            }

            # perform clustering and saves results
            for(scaling in c(TRUE, FALSE)){
                cluster_ranges(dat = df, log_scale = scaling)
            }

            # find cut-off values for each range
            cutoff.ranges <- df %>% 
                dplyr::group_by(Cluster) %>% 
                dplyr::summarise(
                    lower.end = min(Value, na.rm = TRUE),
                    upper.end = max(Value, na.rm = TRUE)
                )

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
                dplyr::arrange(exp, Value)

            ranges.df <- ranges.df %>% 
                dplyr::group_by(exp, Cluster) %>% 
                dplyr::mutate(
                    Cluster = ifelse(Value == min(Value), 
                        paste(Cluster, "-1", sep = ""), 
                        paste(Cluster, "-2", sep = ""))) %>% 
                tidyr::spread(Cluster, Value) %>% 
                dplyr::select(
                    stringr::str_sort(names(.), numeric = TRUE)
                )

            fwrite(
                x = ranges.df, 
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
                        category_sub, 
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
            df.dist <- dist(1-df.hc) %>% suppressWarnings()
            df.dendro <- as.dendrogram(hclust(df.dist, method = "complete"))
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
                df.dendro <- as.dendrogram(hclust(df.dist, method = "complete"))
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
