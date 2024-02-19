FitCurves <- R6::R6Class(
    classname = "FitCurves",
    public = list(
        #' @field bp_dir Numeric vector to fix the seed.
        seed = 1234,

        #' @field chr Numeric vector of chromosome number.
        chr = NULL,

        #' @field which_exp_ind Numeric vector of the index of the experiment 
        #' to process from org_file. If NULL, all exp will be processed.
        which_exp_ind = NULL,

        #' @field results List of results stored if return_vals is TRUE.
        results = NULL,

        initialize = function(seed, chr, which_exp_ind, return_vals, rmsd_values){
            if(!missing(seed)) self$seed <- seed
            if(!missing(chr)) self$chr <- chr
            if(!missing(which_exp_ind)) self$which_exp_ind <- which_exp_ind
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
                    private$to_cluster_curves <- TRUE
                }
                private$len_of_loop <- self$which_exp_ind
            } 
        },

        #' @description
        #' Generate plots from RMSD calculations for each experiment.
        #' @param k size of kmer to plot and fit Gaussian curves to.
        #' @param annot_plots if TRUE, will add 95% C.I. of each curve onto plot.
        #' @param per_chr if TRUE, will fit curves to each chromosome. 
        #'  Default is to take average RMSD per position across chromosomes.
        #' @return None.
        generate_rmsd_plots = function(k = c(2,4,6,8), annot_plots = TRUE, per_chr = FALSE){
            start.time <- Sys.time()
            cur.msg <- "Generating RMSD plots and fitting Gaussian curves"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")

            private$per_chr <- per_chr

            if(is.null(private$rmsd_values)){
                for(i in private$len_of_loop){
                    private$bp_exp <- paste0(
                        private$org_file[i, `Fragmentation type`], "/",
                        private$org_file[i, `Experiment folder`]
                    )
                    private$category <- private$org_file[i, Category_Main]
                    private$category_general <- private$org_file[i, Category_general]

                    cur.msg <- paste0("Processing ", match(i, private$len_of_loop), "/", 
                                    length(private$len_of_loop),
                                    ". ", private$bp_exp)
                    cat(cur.msg, "\n", sep = "")

                    # if k is not specified, perform calculations on c(2,4,6,8)
                    if(length(k) == 1){
                        if(k == 8) private$fits <- c("kmer_8" = 3)
                        if(k == 6) private$fits <- c("kmer_6" = 3)
                        if(k == 4) private$fits <- c("kmer_4" = 3)
                        if(k == 2) private$fits <- c("kmer_2" = 3)
                    } else {
                        private$fits <- c(
                            "kmer_2" = 3, 
                            "kmer_4" = 3, 
                            "kmer_6" = 3, 
                            "kmer_8" = 3
                        )
                    }

                    if(private$per_chr){
                        for(chr in self$chr){
                            t1 <- Sys.time()
                            cur.msg <- paste0("Processing for chromosome ", chr)
                            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                            cat(cur.msg, l, sep = "")

                            private$plot_rmsd(
                                annot_plots = annot_plots,
                                chr = chr
                            )

                            total.time <- Sys.time() - t1
                            cat("DONE! --", signif(total.time[[1]], 2), 
                                attr(total.time, "units"), "\n")
                        }
                    } else {
                        t1 <- Sys.time()
                        cur.msg <- "Processing for all chromosomes"
                        l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                        cat(cur.msg, l, sep = "")

                        private$plot_rmsd(annot_plots = annot_plots)

                        total.time <- Sys.time() - t1
                        cat("DONE! --", signif(total.time[[1]], 2), 
                            attr(total.time, "units"), "\n")
                    }
                }
            } else {
                if(k == 8) private$fits <- c("kmer_8" = 3)
                if(k == 6) private$fits <- c("kmer_6" = 3)
                if(k == 4) private$fits <- c("kmer_4" = 3)
                if(k == 2) private$fits <- c("kmer_2" = 3)

                if(private$per_chr){
                    for(chr in self$chr){
                        t1 <- Sys.time()
                        cur.msg <- paste0("Processing for chromosome ", chr)
                        l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                        cat(cur.msg, l, sep = "")

                        private$plot_rmsd(
                            annot_plots = annot_plots,
                            chr = chr
                        )

                        total.time <- Sys.time() - t1
                        cat("DONE! --", signif(total.time[[1]], 2), 
                            attr(total.time, "units"), "\n")
                    }
                } else {
                    t1 <- Sys.time()
                    cur.msg <- "Processing for all chromosomes"
                    l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                    cat(cur.msg, l, sep = "")

                    private$plot_rmsd(annot_plots = annot_plots)

                    total.time <- Sys.time() - t1
                    cat("DONE! --", signif(total.time[[1]], 2), 
                        attr(total.time, "units"), "\n")
                }
            }

            if(private$to_cluster_curves){
                # cluster rmsd tracts into ranges
                for(kmer in k){                   
                    private$cluster_curves_into_ranges(kmer)
                }
            }

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

        #' @field return_vals Boolean. If TRUE, will return values instead of saving as file.
        return_vals = FALSE,

        #' @field rmsd_values Numeric vector of RMSD values. 
        rmsd_values = NULL,

        #' @field per_chr if TRUE, will fit curves to each chromosome. 
        #'  Default is to take average RMSD per position across chromosomes.
        per_chr = FALSE,

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
        #' @param annot_plots if TRUE, will add 95% C.I. of each curve onto plot.
        #' @param chr Numeric vector of chromosome number.
        #' @return None.
        plot_rmsd = function(annot_plots, chr = 1:22){
            if(is.null(private$rmsd_values)){
                # load data sets
                str_pat <- ifelse(
                    private$per_chr,
                    paste0("^rmsd_", names(private$fits)),
                    paste0(
                        "^chr(?:[1-9]|1[0-9]|2[0-2])_rmsd_",
                        names(private$fits)
                    )
                )
                files <- list.files(
                    path = paste0("../data/", private$bp_exp),
                    pattern = str_pat,
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

                    if(!private$per_chr){
                        out <- as_tibble(out) %>% 
                            dplyr::mutate(
                                chr = paste0("chr", x)
                            )
                    }

                    out <- as_tibble(out) %>% 
                        dplyr::mutate(
                            kmer = as.factor(file.name[x]), 
                            x = -limits:(nrow(out)-limits-1)
                        ) %>% 
                        dplyr::rename(y = value)
                    return(out)
                })
                data.sets <- do.call(rbind, data.sets)

                if(!private$per_chr){
                    data.sets <- data.sets %>% 
                        dplyr::group_by(x) %>% 
                        dplyr::summarise(y = mean(y, na.rm = TRUE))

                    limits <- nrow(data.sets)/2-1
                    data.sets <- data.sets %>% 
                        dplyr::mutate(
                            kmer = as.factor(file.name[1]), 
                            x = -limits:(nrow(data.sets)-limits-1)
                        )
                }

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
                        theme_classic() + 
                        theme(text = element_text(size = 15)) + 
                        labs(
                            x = "Position away from breakpoint, bp",
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
                    if(!private$per_chr){
                        chr <- "all"
                        height <- 5
                        width <- 6
                    }
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
                if(private$per_chr){
                    limits <- length(private$rmsd_values[[paste0("chr", chr)]])/2-1
                    data.sets <- tibble(
                        x = -limits:(length(private$rmsd_values[[paste0("chr", chr)]])-limits-1),
                        y = private$rmsd_values[[paste0("chr", chr)]],
                        kmer = names(private$fits)
                    )
                } else {
                    # for(i in 1:length(files)){
                    #     print(i)
                    #     rmsd_values <- as.data.table(data.sets)
                    #     temp <- rmsd_values[chr == paste0("chr", i)]
                    #     private$rmsd_values[[paste0("chr", i)]] <- temp[["y"]]
                    # }

                    data.sets <- lapply(1:length(names(private$rmsd_values)), function(chr){
                        out <- private$rmsd_values[[names(private$rmsd_values)[chr]]]
                        limits <- length(out)/2-1

                        out <- as_tibble(out) %>% 
                            dplyr::mutate(
                                chr = paste0("chr", chr),
                                kmer = as.factor(names(private$fits[1])), 
                                x = -limits:(length(out)-limits-1)
                            ) %>% 
                            dplyr::rename(y = value)

                        return(out)
                    })
                    data.sets <- do.call(rbind, data.sets)

                    data.sets <- data.sets %>% 
                        dplyr::group_by(x) %>% 
                        dplyr::summarise(y = mean(y, na.rm = TRUE))

                    limits <- nrow(data.sets)/2-1
                    data.sets <- data.sets %>% 
                        dplyr::mutate(
                            kmer = as.factor(names(private$fits)[1]), 
                            x = -limits:(nrow(data.sets)-limits-1)
                        )
                }
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

                    output <- private$make_plot(
                        dat = dat, 
                        k = names(private$fits[k]), 
                        curve.vals = curvefits, 
                        annot_plots = annot_plots
                    )
                    p[[k]] <- output$fit_plot
                    output$fit_plot <- NULL

                    # percent contribution of each gaussian curve towards breakage
                    curve_names <- c("curve.one" = 1, "curve.two" = 2, "curve.three" = 3)
                    df <- as_tibble(t(output$curve_contributions)*100)
                    colnames(df) <- names(curve_names[1:ncol(df)])
                    if(ncol(df) == 1) df <- dplyr::mutate(df, curve.two = NA_real_)
                    if(ncol(df) == 2) df <- dplyr::mutate(df, curve.three = NA_real_)

                    # range of influence based on 95 percent confidence intervals
                    output.CIlst <- output$conf_int
                    curves <- output.CIlst[seq(2, length(output.CIlst), 2)]-
                              output.CIlst[seq(1, length(output.CIlst), 2)]
                    curves <- c(curves, rep(NA_real_, 3-length(curves)))

                    # sd of curves
                    SD <- as.numeric(output$curve_sds)
                    SD <- c(SD, rep(NA_real_, 3-length(SD)))

                    # contribution of the peak intensity of each gaussian curves
                    peak.intensity <- as.numeric(output$peak_contributions)*100
                    peak.intensity <- c(peak.intensity, rep(NA_real_, 3-length(peak.intensity)))

                    # combine results
                    df <- rbind(df, curves, SD, peak.intensity)
                    df <- df %>% 
                        dplyr::mutate(
                            exp = private$bp_exp,
                            category = private$category,
                            category_general = private$category_general,
                            rowid = c("contribution", "ranges", "SD", "peak.intensity"), 
                            .before = 1
                        )

                    # save the individual points of the gaussian curve
                    curve.vals <- lapply(1:unname(output$curves), function(x){
                        temp <- tibble(vals = output$model_fit[[x]]) %>% 
                            dplyr::mutate(
                                col.label = dplyr::case_when(
                                    x == 1 ~ "curve.one",
                                    x == 2 ~ "curve.two",
                                    x == 3 ~ "curve.three"
                                )
                            ) %>% 
                            as.data.table()
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
                        df <- df[, 1:(4+curvefits$curves)]
                    }
                    all.gc <- as.data.table(output$gaussian_linear_comb)

                    if(private$return_vals | private$per_chr){                        
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
            }

            if(is.null(private$rmsd_values)){
                height <- ifelse(length(files) == 1, 5, 12)
                width <- ifelse(length(files) == 1, 7, 5)
                if(!private$return_vals){
                    if(!private$per_chr){
                        chr <- "all"
                        height <- 5
                        width <- 6
                    }
                    pdf(
                        file = paste0("../figures/", private$bp_exp, 
                                    "/chr", chr, "_RMSD_all_GMMFit.pdf"),
                        height = height, width = width
                    )
                    p1 <- do.call(gridExtra::grid.arrange, c(p, ncol = 1))
                    pic.saved <- dev.off()
                }
            }
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

            # fit each Gaussian curve separately to data
            draw.from.gaussian <- function(xs, C, SD){
                return(C*exp(-(xs-0)^2/(2*SD^2))+min(dat$y, na.rm = TRUE))
            }

            calc_contributions <- function(dat, summary.fit.params){
                num_curves <- length(summary.fit.params)/2
                curve_labels <- paste0("C", 1:num_curves)
                sd_labels <- paste0("sd", 1:num_curves)
                
                # Draw Gaussian curves and calculate peak intensity for each curve
                fits <- lapply(1:num_curves, function(i){
                    draw.from.gaussian(
                        xs = dat$x,
                        C = summary.fit.params[curve_labels[i]],
                        SD = summary.fit.params[sd_labels[i]]
                    )
                })
                if(length(fits) > 1){
                    peak.intensities <- sapply(2:length(fits), function(fit){
                        max(fits[[fit]], na.rm = TRUE)-min(dat$y, na.rm = TRUE)
                    })
                    peak.intensities <- c(max(fits[[1]], na.rm = TRUE), peak.intensities)
                } else {
                    peak.intensities <- max(fits[[1]], na.rm = TRUE)
                }
                
                # Calculate integral (area under the curve) for each Gaussian
                integral.curves <- summary.fit.params[curve_labels] * 
                                    summary.fit.params[sd_labels] * 
                                    sqrt(2*pi)

                # Calculate total area under all curves
                integral.all.curves <- sum(integral.curves)
                
                # Calculate percentage contribution of each curve
                curve.contributions <- integral.curves/integral.all.curves

                # Calculate peak intensity contributions
                peak.intensity.all <- sum(peak.intensities)
                peak.contributions <- peak.intensities/peak.intensity.all
                
                # Check if any peak contributions are less than 5%
                return(
                    list(
                        model_fit = fits,
                        peak_contributions = peak.contributions,
                        curve_contributions = curve.contributions,
                        curve_sds = summary.fit.params[sd_labels],
                        error = any(peak.contributions*100 < 5)
                    )
                )
            }

            fit_gaussians <- function(dat, curves){
                # Create a dynamic grid search base that adheres to the condition sd1 < sd2 < sd3
                set.seed(self$seed)
                num_samples <- 1000
                upper_max <- abs(max(x))-1

                sample_space <- seq(from = 0.1, to = ceiling(upper_max/4), by = 0.01)
                sd_ranges <- lapply(1:curves, function(i){
                    # # Create weights that bias towards smaller values for fewer curves
                    # weights <- (1:length(sample_space))^(-log(i+1))
                    sample(sample_space, num_samples, replace = FALSE)
                })
                grid_search <- do.call(cbind, sd_ranges)
                grid_search <- t(apply(grid_search, 1, sort))
                grid_search <- grid_search[order(grid_search[,1]),]
                grid_search <- as.data.frame(grid_search)
                colnames(grid_search) <- paste0("sd", 1:curves)
                grid_search <- grid_search[, 1:curves]

                # # Create a dynamic grid search base that adheres to the condition sd1 < sd2 < sd3
                # sd_ranges <- lapply(1:curves, function(i){
                #     seq(from = i*0.2, to = 150, length.out = 500)
                # })
                # grid_search <- as.data.frame(do.call(cbind, sd_ranges))
                # colnames(grid_search) <- paste0("sd", 1:curves)
                # grid_search <- grid_search[, 1:curves]
                best_fit <- list(fit = NULL, rss = Inf, curves = curves)

                repeat {
                    # Build the formula string based on the number of curves
                    formula_parts <- paste0("C", 1:curves, " * exp(-(x-0)^2/(2 * sd", 1:curves, "^2))")
                    formula_str <- paste("y ~", paste(formula_parts, collapse=" + "), "+", "min(y)")
                    if(is.null(dim(grid_search))){
                        grid_search <- data.frame("sd1" = grid_search)
                    }
                    grid_search$rss <- rep(NA, nrow(grid_search))
                    sd_names <- paste0("sd", 1:curves)

                    # Attempt to fit the model using grid search
                    for(index in 1:nrow(grid_search)){
                        coefficients_start <- setNames(rep(10, curves), paste0("C", 1:curves))        
                        sds_start <- as.numeric(grid_search[index, sd_names])
                        start_list <- c(coefficients_start, sds_start)
                        names(start_list) <- c(paste0("C", 1:curves), sd_names)

                        tryCatch({
                            fit <- nls(
                                formula_str, 
                                data = dat, 
                                start = start_list, 
                                control = nls.control(
                                    maxiter = 50000, 
                                    tol = 1e-05, 
                                    warnOnly = TRUE
                                ),
                                lower = rep(c(0, 0), each = curves),
                                upper = rep(c(NULL, upper_max), each = curves),
                                trace = FALSE,
                                algorithm = "port"
                            ) %>% suppressWarnings()

                            # Check fit parameters
                            if(class(fit) == "nls"){
                                check.params <- fit$m$getPars()
                                # check for zero coefficients
                                any.zero.coefs <- any(check.params[grepl(
                                    "^C", x = names(check.params)
                                )] == 0)
                                # check for any out-of-bound sds 
                                any.oob.sds <- any(check.params[grepl(
                                    "^sd", x = names(check.params)
                                )] <= 1)

                                # check for peak contribution threshold
                                check.contr <- calc_contributions(
                                    dat = dat, 
                                    summary.fit.params = check.params
                                )

                                # print(paste0(
                                #     "Curves:       ", curves, ". ",
                                #     "Index:        ", index,  ". ",
                                #     "Zero Cs:      ", any.zero.coefs,  ". ",
                                #     "OOB SDs:      ", any.oob.sds,  ". ",
                                #     "Curve contr:  ", check.contr$error,  ". "
                                # ))

                                if(any(c(any.zero.coefs, any.oob.sds, check.contr$error))){
                                    next
                                }
                            }

                            # If fit is successful and parameters are valid, calculate RSS
                            grid_search$rss[index] <- sum(residuals(fit)^2)
                            if(grid_search$rss[index] < best_fit$rss){
                                best_fit <- list(
                                    fit = fit, 
                                    rss = grid_search$rss[index], 
                                    curves = curves
                                )
                            }

                            # early stopping criterion if rss doesn't improve
                            last_rss <- tail(grid_search$rss[!is.na(grid_search$rss)], n = 5)
                            if(length(last_rss) == 5){
                                last_rss <- signif(last_rss, digits = 8)
                                cur_rss <- signif(grid_search$rss[index], digits = 8)
                                are_identical <- length(unique(last_rss)) == 1
                                if(are_identical & identical(unique(last_rss), cur_rss)){
                                    break
                                }
                            }
                        }, error = function(e){
                            NULL
                        })
                    }

                    # if a fit is found, break the loop
                    if(!is.null(best_fit$fit)) break

                    # If no fit is found, decrease the number of curves and try again
                    curves <- curves-1
                    if(curves < 1){
                        break
                    } else {
                        # Update grid_search to match the new number of curves
                        grid_search <- grid_search[, paste0("sd", 1:curves), drop = FALSE]
                    }
                }
                return(best_fit)
            }

            fit_results <- fit_results_3 <- fit_gaussians(
                dat = dat, curves = nr.of.curves
            )
            if(fit_results_3$curves == nr.of.curves){
                fit_results_2 <- fit_gaussians(
                    dat = dat, curves = nr.of.curves-1
                )
                if(fit_results_2$rss < fit_results_3$rss){
                    fit_results <- fit_results_2
                }
            }
            model <- fit_results$fit
            rss <- fit_results$rss
            curves <- fit_results$curves

            init.list <- vector(mode = "list", length = 4)
            init.list[[4]] <- model

            # extract each parameter
            summary.fit.params <- model$m$getAllPars()

            # get contributions
            model.fits <- calc_contributions(
                dat = dat, 
                summary.fit.params = summary.fit.params
            )
            model.fits$model <- model
            model.fits$rss <- rss
            model.fits$curves <- curves
            model.fits$coefs <- summary.fit.params

            return(model.fits)
        },

        #' @description
        #' Make the plot of RMSD values with curves fitted to those values.
        #' @param dat Tibble of the RMSD values per k-mer.
        #' @param k Numeric vector of the k-mer.
        #' @param curve.vals Numeric vector of the values for the fitted curve to plot.
        #' @param annot_plots if TRUE, will add 95% C.I. of each curve onto plot.
        #' @return list of plots.
        make_plot = function(dat, k, curve.vals, annot_plots){
            # extract standard deviation values to compute 95% confidence intervals
            coefs <- curve.vals$coefs
            st.devs <- unname(coefs[grepl("sd", names(coefs))])
            CI.lst <- numeric(length = curve.vals$curves*2)
            CI.lst[1] <- -1.96*st.devs[[1]]
            CI.lst[2] <- 1.96*st.devs[[1]]
            limits <- 500

            if(curve.vals$curves > 1){
                CI.lst[3] <- -1.96*st.devs[[2]]
                CI.lst[4] <- 1.96*st.devs[[2]]
            } 

            if(curve.vals$curves == 3){
                CI.lst[5] <- -1.96*st.devs[[3]]
                CI.lst[6] <- 1.96*st.devs[[3]]
            }
            y.pos <- max(dat$y)

            # get data into the right format
            rmsd_vals <- lapply(1:curve.vals$curves, function(x){
                rmsd_vals <- curve.vals$model_fit[[x]]
                limits <- length(rmsd_vals)/2-1
                x_vals <- -limits:(length(rmsd_vals)-limits-1)

                rmsd_vals <- tibble(y = rmsd_vals) %>% 
                    dplyr::mutate(
                        curve = dplyr::case_when(
                            x == 1 ~ "curve.one",
                            x == 2 ~ "curve.two",
                            x == 3 ~ "curve.three"
                        ),
                        lower_CI = dplyr::case_when(
                            x == 1 ~ CI.lst[1],
                            x == 2 ~ CI.lst[3],
                            x == 3 ~ CI.lst[5]
                        ),
                        upper_CI = dplyr::case_when(
                            x == 1 ~ CI.lst[2],
                            x == 2 ~ CI.lst[4],
                            x == 3 ~ CI.lst[6]
                        ),
                        kmer = names(curve.vals$curves),
                        x = x_vals
                    )
                return(rmsd_vals)
            })
            rmsd_vals <- do.call(rbind, rmsd_vals)
            rmsd_vals <- rmsd_vals %>% 
                dplyr::mutate(
                    annot_CI = dplyr::case_when(
                        curve.vals$curves == 1 ~ paste0(
                            "95% CIs:\n",
                            "R: [", signif(CI.lst[1], 3), ",", 
                            signif(CI.lst[2], 3), "]"
                        ),
                        curve.vals$curves == 2 ~ paste0(
                            "95% CIs:\n",
                            "R: [", signif(CI.lst[1], 3), ",", 
                            signif(CI.lst[2], 3), "]\n",
                            "B: [", signif(CI.lst[3], 3), ",", 
                            signif(CI.lst[4], 3), "]"
                        ),
                        curve.vals$curves == 3 ~ paste0(
                            "95% CIs:\n",
                            "R: [", signif(CI.lst[1], 3), ",", 
                            signif(CI.lst[2], 3), "]\n",
                            "B: [", signif(CI.lst[3], 3), ",", 
                            signif(CI.lst[4], 3), "]\n",
                            "G: [", signif(CI.lst[5], 3), ",", 
                            signif(CI.lst[6], 3), "]"
                        )
                    )
                ) %>% 
                dplyr::select(
                    x, y, kmer, curve, 
                    lower_CI, upper_CI, annot_CI
                )

            gaussian_linear_comb <- rmsd_vals %>% 
                dplyr::group_by(x) %>% 
                dplyr::mutate(
                    y = sum(y),
                    y = dplyr::case_when(
                        curve.vals$curves == 3 ~ y-2*min(dat$y),
                        curve.vals$curves == 2 ~ y-min(dat$y),
                        curve.vals$curves == 1 ~ y
                    )
                ) %>% 
                dplyr::filter(curve == "curve.one") %>% 
                dplyr::mutate(
                    curve = "curve.all",
                    hex = "#f8b195"
                )

            dat <- dat %>% 
                dplyr::mutate(
                    curve = "RMSD",
                    lower_CI = NA_integer_,
                    upper_CI = NA_integer_,
                    annot_CI = NA_character_
                ) %>% 
                dplyr::select(
                    x, y, kmer, curve, 
                    lower_CI, upper_CI, annot_CI
                ) %>% 
                dplyr::bind_rows(rmsd_vals, gaussian_linear_comb) %>% 
                dplyr::mutate(
                    hex = dplyr::case_when(
                        curve == "RMSD" ~ "#000000",
                        curve == "curve.one" ~ "#619B61",
                        curve == "curve.two" ~ "#355c7d",
                        curve == "curve.three" ~ "#BB357E",
                        curve == "curve.all" ~ "#f8b195"
                    ),
                    alpha = dplyr::case_when(
                        curve == "RMSD" ~ 0.3,
                        TRUE ~ 1
                    ),
                    lw = dplyr::case_when(
                        curve == "RMSD" ~ 1,
                        curve == "curve.one" ~ 1.2,
                        curve == "curve.two" ~ 1.2,
                        curve == "curve.three" ~ 1.2,
                        curve == "curve.all" ~ 2
                    ),
                    curve = factor(curve, levels = c(
                            "RMSD",
                            "curve.one", 
                            "curve.two", 
                            "curve.three", 
                            "curve.all"
                        )
                    ),
                    kmer = stringr::str_replace_all(
                        string = stringr::str_to_title(kmer),
                        pattern = "_",
                        replacement = " "
                    )
                )

            # Plot the data with the model superimposed
            fit.plot <- dat %>%
                ggplot(aes(x = x, y = y, col = hex)) + 
                geom_line(aes(alpha = alpha), linewidth = dat$lw) +
                facet_wrap(vars(kmer)) + 
                scale_color_identity() + 
                theme_bw() + 
                theme_classic() + 
                theme(
                    text = element_text(size = 15),
                    legend.position = "none"
                ) + 
                labs(
                    x = "Position away from breakpoint, bp",
                    y = "RMSD"
                )

            if(annot_plots){
                filtered_dat_lower <- dat %>% 
                    dplyr::filter(!is.na(lower_CI))
                filtered_dat_upper <- dat %>% 
                    dplyr::filter(!is.na(upper_CI))
                annotation_text <- unique(dat$annot_CI[!is.na(dat$annot_CI)])
                
                fit.plot <- fit.plot + 
                    geom_vline(
                        data = filtered_dat_lower,
                        aes(xintercept = lower_CI, color = hex), 
                        linetype = "dashed",
                        alpha = 0.7
                    ) +
                    geom_vline(
                        data = filtered_dat_upper,
                        aes(xintercept = upper_CI, color = hex), 
                        linetype = "dashed",
                        alpha = 0.7
                    ) +
                    geom_text(
                        aes(
                            x = -Inf, y = Inf,
                            label = annotation_text, 
                            hjust = -0.1, vjust = 1.1
                        ),
                        size = 4,
                        color = "black"
                    )
            }
            curve.vals$gaussian_linear_comb <- gaussian_linear_comb$y
            curve.vals$fit_plot <- fit.plot
            curve.vals$conf_int <- CI.lst
            return(curve.vals)
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
                        bins = 100
                    ) +
                    theme_bw() + 
                    theme_classic() + 
                    theme(text = element_text(size = 15)) + 
                    scale_fill_identity() + 
                    labs(
                        subtitle = "Clustering based on true values",
                        x = "", y = ""
                    )

                if(log_scale){
                    dat.plot <- dat.plot + 
                        labs(
                            subtitle = "Clustering based on log2-scaled values",
                            x = "Sequence context range, bp",
                        )
                }
                    
                dir.create(
                    path = "../figures/ranges",
                    showWarnings = FALSE,
                    recursive = TRUE
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

                return(dat.plot)
            }

            # perform clustering and saves results
            clustering_plots <- lapply(c(TRUE, FALSE), function(scaling){
                return(cluster_ranges(dat = df, log_scale = scaling))
            })

            pdf(
                file= paste0(
                    "../figures/ranges/kmer_", kmer, 
                    "_LOG-and-raw_flattened_ranges_histogram.pdf"
                ),
                width = 7, height = 10
            )
            gridExtra::grid.arrange(
                clustering_plots[[2]], 
                clustering_plots[[1]], 
                ncol = 1
            )
            save.plot <- dev.off()

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
                tidyr::pivot_wider(
                    names_from = Cluster,
                    values_from = Value
                ) %>% 
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