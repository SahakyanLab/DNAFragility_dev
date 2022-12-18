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
                private$len_of_loop <- 1:nrow(private$org_file)
            } else {
                private$len_of_loop <- self$which_exp_ind
            }
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
                t1 <- Sys.time()
                cur.msg <- paste0("Generating RMSD plots for exp ", 
                                  i, " of ", nrow(private$org_file))
                l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                cat(paste0(cur.msg, l))

                private$bp_exp <- paste0(
                    private$org_file[i, `Fragmentation type`], "/",
                    private$org_file[i, `Experiment folder`]
                )
                private$category <- private$org_file[i, Category_Main]
                private$category_sub <- private$org_file[i, Category_sub]
                private$fits <- c("kmer_4" = 3, "kmer_6" = 3, "kmer_8" = 3)
                private$plot_rmsd()

                total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")
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
            cat("Final time taken:", signif(final.t[[1]], digits = 3), 
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
            df <- fread(
                "../../data/org_file.csv",
                showProgress = FALSE
            )
            private$org_file <- df[`DSB Map` == "TRUE"]
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
            file.name <- paste0("kmer_", seq(from = 2, to = 8, by = 2))
            files <- stringr::str_sort(files, numeric = TRUE)

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

            if(private$auto_fit){
                # Overall RMSD plots
                plots <- data.sets %>%
                    ggplot(aes(x = x, y = y)) + 
                        geom_line(size = 0.8) + 
                        facet_wrap(~kmer, ncol = 2, scales = "free_y") + 
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
                    height = 7, 
                    width = 7
                )
            }

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
                        file = paste0("../data/", private$bp_exp, 
                                      ifelse(private$auto_fit, "/key", "/new"),
                                      "_stats_", names(private$fits[k]), ".csv"),
                        row.names = FALSE
                    )
                }
                # if(k == 1) print(p) # temp. avoid an unsolved ggplot bug
            }

            ggsave(
                filename = paste0("../figures/", private$bp_exp, 
                                  "/chr", self$chr, "_RMSD_all_",
                                  ifelse(private$auto_fit, 
                                  "GMMFit.pdf", "GMMFit_optimised.pdf")),
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

            # try_curvefit <- function(C.value, nr.of.curves){
            #     # model-based clustering to obtain starting st.dev values
            #     # for curve fittings
            #     set.seed(1234)
            #     best_fit <- mclust::Mclust(
            #         dat.norm, 
            #         G = nr.of.curves, 
            #         verbose = FALSE
            #     )
            #     params <- matrix(
            #         best_fit$parameters$variance$sigma, 
            #         ncol = nr.of.curves
            #     )
            #     params.vars <- ceiling(sqrt(params[1,])/2)
            #     params.sds <- sort(params.vars)

            #     if(nr.of.curves == 2){
            #         fit_once <- nls(
            #             y ~ (C1 * exp(-(x-0)**2/(2 * sigma1**2)) +
            #                  C2 * exp(-(x-0)**2/(2 * sigma2**2)) +
            #                  min(y)),
            #             data=dat,
            #             start=list(C1=C.value, sigma1=params.sds[1],
            #                        C2=C.value, sigma2=params.sds[2]
            #                        ),
            #             control = nls.control(
            #                 maxiter=50000, 
            #                 tol=1e-05, 
            #                 warnOnly = TRUE
            #             ),
            #             lower=rep(c(0, 0), 2),
            #             upper=rep(c(NULL, 150), 2),
            #             trace=FALSE,
            #             algorithm="port") %>% 
            #             suppressWarnings()
            #     } else if(nr.of.curves == 3){
            #         fit_once <- nls(
            #             y ~ (C1 * exp(-(x-0)**2/(2 * sigma1**2)) +
            #                  C2 * exp(-(x-0)**2/(2 * sigma2**2)) +
            #                  C3 * exp(-(x-0)**2/(2 * sigma3**2)) +
            #                 min(y)),
            #             data=dat,
            #             start=list(C1=C.value, sigma1=params.sds[1],
            #                       C2=C.value, sigma2=params.sds[2],
            #                       C3=C.value, sigma3=params.sds[3]
            #                       ),
            #             control = nls.control(
            #                 maxiter=50000, 
            #                 tol=1e-05, 
            #                 warnOnly = TRUE
            #             ),
            #             lower=rep(c(0, 0), 3),
            #             upper=rep(c(NULL, 150), 3),
            #             trace=FALSE,
            #             algorithm="port") %>% 
            #             suppressWarnings()
            #     }
            #     return(fit_once)
            # }

            # fit <- tryCatch({
            #     try_curvefit(
            #         C.value = C.value, 
            #         nr.of.curves = nr.of.curves
            #     )
            # }, error = function(e) return(2))
            # if(class(fit) != "nls"){
            #     nr.of.curves <<- nr.of.curves-1
            #     fit <- try_curvefit(
            #         C.value = C.value, 
            #         nr.of.curves = nr.of.curves
            #     )
            # }

            grid_search <- data.frame(
                sd1 = seq(2, 18, 1),
                sd2 = seq(8, 40, 2),
                sd3 = seq(20, 70, 3)
            )

            if(nr.of.curves == 2){
                try_fit <- function(C.value = 10, sd1, sd2){
                    nls(
                    y ~ (C1 * exp(-(x-0)**2/(2 * sigma1**2)) +
                            C2 * exp(-(x-0)**2/(2 * sigma2**2)) +
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
                    upper=rep(c(NULL, 150), 2),
                    trace=FALSE,
                    algorithm="port") %>% 
                    suppressWarnings()
                }

                index <- 1
                while(index <= nrow(grid_search)){
                    fit <- tryCatch({
                        try_fit(
                            C.value = 10,
                            sd1 = grid_search$sd1[index],
                            sd2 = grid_search$sd2[index]
                        )
                    }, error = function(e) return(2))
                    if(class(fit) == "nls"){
                        check.params <- fit$m$getAllPars()
                        # check for zero coefficients
                        any.zero.coefs <- any(check.params[grepl(
                            pattern = "C", x = names(check.params)
                        )] == 0)
                        # check for sigma too small
                        any.small.sigmas <- any(check.params[grepl(
                            pattern = "sigma", x = names(check.params)
                        )] < (grid_search[1,"sd1"]/2))
                        if(any.zero.coefs | any.small.sigmas){
                            index <- index+1
                        } else {
                            break
                        }
                    } else {
                        index <- index+1
                    }
                }
            } else if(nr.of.curves == 3){
                try_fit <- function(C.value = 10, sd1, sd2, sd3){
                    nls(
                    y ~ (C1 * exp(-(x-0)**2/(2 * sigma1**2)) +
                            C2 * exp(-(x-0)**2/(2 * sigma2**2)) +
                            C3 * exp(-(x-0)**2/(2 * sigma3**2)) +
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
                    upper=rep(c(NULL, 150), 3),
                    trace=FALSE,
                    algorithm="port") %>% 
                    suppressWarnings()
                }

                index <- 1
                while(index <= nrow(grid_search)){
                    fit <- tryCatch({
                        try_fit(
                            C.value = 10,
                            sd1 = grid_search$sd1[index],
                            sd2 = grid_search$sd2[index],
                            sd3 = grid_search$sd3[index]
                        )
                    }, error = function(e) return(2))
                    if(class(fit) == "nls"){
                        check.params <- fit$m$getAllPars()
                        # check for zero coefficients
                        any.zero.coefs <- any(check.params[grepl(
                            pattern = "C", x = names(check.params)
                        )] == 0)
                        # check for sigma too small
                        any.small.sigmas <- any(check.params[grepl(
                            pattern = "sigma", x = names(check.params)
                        )] < (grid_search[1,"sd1"]/2))
                        if(any.zero.coefs | any.small.sigmas){
                            index <- index+1
                        } else {
                            break
                        }
                    } else {
                        index <- index+1
                    }
                }
            }

            # error capture
            init.list <- vector(mode = "list", length = 4)
            init.list[[4]] <- fit
            if(class(fit) != "nls") return(init.list) 

            # extract each parameter
            summary.fit.params <- fit$m$getAllPars()

            # check for zero coefficients
            any.zero.coefs <- any(summary.fit.params[grepl(
                pattern = "C", x = names(summary.fit.params)
            )] == 0)
            if(any.zero.coefs) return(init.list[[4]] <- 2)

            # check for sigma too small
            any.small.sigmas <- any(summary.fit.params[grepl(
                pattern = "sigma", x = names(summary.fit.params)
            )] < grid_search[1, "sd1"])
            if(any.zero.coefs) return(init.list[[4]] <- 2)
            
            # fit each Gaussian curve separately to data
            draw.from.gaussian <- function(xs, C, SD, miny) {
                return(C * exp(-(xs-0)**2/(2 * SD**2)) + min(miny))
            }

            # integral functions for curve 1
            curve.one <- function(x, C, SD){
                return(
                    (C*exp(-(x-0)**2/(2*SD**2)))
                )
            }

            # integral functions for curve 2
            curve.two <- function(x, C, SD){
                return(
                    (C*exp(-(x-0)**2/(2*SD**2)))
                )
            }

            # integral functions for curve 3
            curve.three <- function(x, C, SD){
                return(
                    (C*exp(-(x-0)**2/(2*SD**2)))
                )
            }

            # integral functions for all curves
            all.curves <- function(x, nr.of.curves, C1, C2, 
                                   C3=NULL, SD1, SD2, SD3=NULL){
                if(nr.of.curves == 2){
                    return(
                        (C1*exp(-(x-0)**2/(2*SD1**2))+
                        C2*exp(-(x-0)**2/(2*SD2**2))
                        )
                    )
                } else if(nr.of.curves == 3){
                    return(
                        (C1*exp(-(x-0)**2/(2*SD1**2))+
                        C2*exp(-(x-0)**2/(2*SD2**2))+
                        C3*exp(-(x-0)**2/(2*SD3**2))
                        )
                    )
                }
            }
            
            # fit curve 1
            fit_1 <- draw.from.gaussian(
                xs = x,
                C = summary.fit.params["C1"],
                SD = summary.fit.params["sigma1"],
                miny = min(y)
            )
            if(max(fit_1, na.rm = TRUE) <= median(y, na.rm = TRUE)){
                return(init.list[[4]] <- 2)
            }

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
                    miny = min(y)
                )
                if(max(fit_2, na.rm = TRUE) <= median(y, na.rm = TRUE)){
                    return(init.list[[4]] <- 2)
                }

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
                        miny = min(y)
                    )
                    if(max(fit_3, na.rm = TRUE) <= median(y, na.rm = TRUE)){
                        return(init.list[[4]] <- 2)
                    }

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

            if(any(CI.lst[2] > 300 | CI.lst[1] < -300)){
                return(2)
            }

            if(nr.of.curves > 1){
                if(any(1.96*st.devs[[2]] > 300 | -1.96*st.devs[[2]] < -300)){
                    return(2)
                }
                CI.lst[3] <- -1.96*st.devs[[2]]
                CI.lst[4] <- 1.96*st.devs[[2]]
            } 

            if(nr.of.curves == 3){
                if(any(1.96*st.devs[[3]] > 300 | -1.96*st.devs[[3]] < -300)){
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
                    size = 1.2) + 
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

            if(nr.of.curves == 1){
                fit.plot <- fit.plot + 
                geom_text(
                    data = data.frame(
                        xpos = -Inf,
                        ypos = Inf,
                        annotateText = paste0(
                            "Red CI:   ", formatC(signif(CI.lst[1], 3), 3), 
                            " / ", formatC(signif(CI.lst[2], 3), 3)
                        ),
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
                    size = 1.2) +
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
                                size = 1) + 
                        geom_text(
                            data = data.frame(
                                xpos = -Inf,
                                ypos = Inf,
                            annotateText = paste0(
                                "Red CI:   ", formatC(signif(CI.lst[1], 3), 3), " / ", 
                                formatC(signif(CI.lst[2], 3), 3), "\n",
                                "Blue CI:  ", formatC(signif(CI.lst[3], 3), 3), " / ", 
                                formatC(signif(CI.lst[4], 3), 3)
                            ),
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
                            size = 1.2) + 
                        geom_line(
                            data = data.frame(
                                x = dat$x, 
                                y = curve.vals[[1]]+curve.vals[[2]]+curve.vals[[3]]-2*min(dat$y)
                            ),
                            stat = "identity",
                            color = "orange",
                            alpha = 0.8,
                            size = 1) + 
                        geom_text(
                            data = data.frame(
                                xpos = -Inf,
                                ypos = Inf,
                                annotateText = paste0(
                                    "Red CI:   ", formatC(signif(CI.lst[1], 3), 3), " / ", 
                                    formatC(signif(CI.lst[2], 3), 3), "\n",
                                    "Green CI: ", formatC(signif(CI.lst[5], 3), 3), " / ", 
                                    formatC(signif(CI.lst[6], 3), 3), "\n",
                                    "Blue CI:  ", formatC(signif(CI.lst[3], 3), 3), " / ", 
                                    formatC(signif(CI.lst[4], 3), 3)
                                ),
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
        }
    )
)