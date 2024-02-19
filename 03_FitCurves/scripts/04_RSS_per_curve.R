# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(cluster)))
suppressPackageStartupMessages(suppressWarnings(library(mclust)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults
pbapply::pboptions(char = "=")

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

seed <- 1234

#' 1. Import data
df <- tibble(
    class = c(
        "Mechanical",
        "Natural decay",
        "Cell free DNA",
        "Endogenous DSBs",
        "Enzymatic"
    ),
    bp_exp = c(
        "00_Breakage/Ultrasonication/Simons_exp_1",
        "00_Breakage/Ancient_DNA/Altai_Neandertal",
        "00_Breakage/cfDNA/Ovarian_cancer",
        "00_Breakage/sBLISS/K562_Top2_mediated_DSBs/ETO",
        "00_Breakage/Enzymatic/EcoRV_HeLa_cells"
    ),
    hex = c(
        "#8B8000",
        "#9CDC9C",
        "#FFA435",
        "#BFABCB",
        "#BB357E"   
    )
)

results <- lapply(1:nrow(df), function(row){
    cat('Processing: ', df$bp_exp[row], '...\n')

    files <- list.files(
        path = paste0("../data/", df$bp_exp[row]),
        pattern = "^chr(?:[1-9]|1[0-9]|2[0-2])_rmsd_kmer_8",
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
                chr = paste0("chr", x)
            )

        out <- as_tibble(out) %>% 
            dplyr::mutate(
                kmer = as.factor(file.name[x]), 
                x = -limits:(nrow(out)-limits-1)
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
            kmer = as.factor(file.name[1]), 
            x = -limits:(nrow(data.sets)-limits-1)
        )

    fit_gaussians <- function(dat, curves, seed, x){
        # Create a dynamic grid search base that adheres to the condition sd1 < sd2 < sd3
        set.seed(seed)
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
        best_fit <- list(fit = NULL, rss = Inf, curves = curves)

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

                    if(any(c(any.zero.coefs, any.oob.sds))){
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
        return(best_fit)
    }

    # Start the fitting process with the desired number of curves
    kmer_size = "kmer_8"
    dat <- data.sets %>% dplyr::filter(kmer == kmer_size)

    #' 3. Fit curves
    x <- dplyr::pull(dat, x)
    y <- dplyr::pull(dat, y)
    dat <- data.frame(x = x, y = y)

    # average normalised behaviour
    dat.norm <- tibble(
        x = x,
        y_fwd = scale(y, center = TRUE, scale = TRUE),
        y_rev = rev(y_fwd)
    ) 
    dat.norm <- data.frame(
        x = x,
        y = rowMeans(dat.norm[2:3], na.rm = TRUE)
    )
    dat.norm <- dat.norm[1:(ceiling(nrow(dat.norm)/2)),]

    fit_results <- pbapply::pblapply(1:7, function(curves){
        fit_gaussians(dat=dat, curves=curves, seed=seed, x=x)
    })

    #' the goal is to represent a measure of 
    #' improvement which we expect to increase.
    #' Thus we will consider using the maximum rss as 
    #' a baseline as that is the worst performer
    #' and we calculate the percentage change 
    #' relative to this baseline.
    df_fit <- tibble(
            class = df$class[row],
            hex = df$hex[row],
            curves = sapply(fit_results, `[[`, 3),
            rss = sapply(fit_results, `[[`, 2)
        ) %>%
        dplyr::mutate(
            norm_rss = (rss-min(rss))/(max(rss))*100,
            perc_chg_exp = 100-norm_rss
        ) %>%
        dplyr::select(-norm_rss) 

    top_row <- df_fit[1,]
    top_row$curves <- 0
    top_row$perc_chg_exp <- 0
    df_fit <- dplyr::bind_rows(top_row, df_fit)
    return(df_fit)
})
results <- do.call(rbind, results)

p1 <- results %>% 
    ggplot(aes(x = curves, y = perc_chg_exp, col = hex)) + 
    geom_line(linewidth = 1.5) +
    geom_point(size = 2) + 
    geom_vline(
        xintercept = 3, 
        linewidth = 1.2,
        colour = "darkgrey", 
        linetype = "dashed"
    ) + 
    theme_bw() + 
    theme_classic() + 
    scale_color_identity() + 
    scale_x_continuous(breaks = unique(results$curves)) +
    theme(text = element_text(size = 20)) + 
    coord_cartesian(ylim = c(0, 100)) + 
    labs(
        title = "Diminishing returns of additional Gaussian curves fitted",
        x = "Number of Gaussian curves fitted",
        y = "Cumulative change in RSS, %"
    )

ggsave(
    filename = "../figures/RSS_Gaussian_curves.pdf",
    plot = p1,
    height = 8, width = 10
)