suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults
pbapply::pboptions(char = "=")

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

# import exact breakage site overlaps
df_exactbreaks <- fread("../data/exact_breaks/overlap_percent.csv")
print(paste("Mean breakage overlaps:", signif(mean(df_exactbreaks$overlaps) * 100, 3)))
# setorder(df_exactbreaks, -overlaps); df_exactbreaks

# import pearson correlation coefficient of k-mer fragility scores
cor.mat <- lapply(c("zscore", "ratio"), function(action){
    res <- pbapply::pblapply(c(2,4,6,8), function(kmer){
        # import datasets
        files <- list.files(
            path = "../data/kmertone",
            pattern = paste0("score_", kmer, "-mers.csv"),
            full.names = TRUE,
            recursive = TRUE
        )
        files <- files[grepl(pattern = "00_Breakage", x = files)]

        exp.labels <- stringr::str_remove_all(
            string = files,
            pattern = paste0("../data/kmertone/|/score_", 
                                kmer, "-mers.csv")
        )
        org.file <- fread("../../data/org_file.csv")
        org.file[, exp := paste0(`Fragmentation type`, "/", `Experiment folder`)]
        breaktype.labels <- org.file$Category_Main[match(exp.labels, org.file$exp)]

        # get correlation density plot of experiments from same type
        df <- tibble(
            file.name = files,
            exp = exp.labels,
            bp.type = breaktype.labels) %>% 
            dplyr::group_by(bp.type) %>% 
            dplyr::mutate(count = dplyr::n()) %>% 
            dplyr::filter(count > 1) %>% 
            dplyr::select(-count) %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate(
                exp = stringr::str_remove_all(
                    string = exp,
                    pattern = "00_Breakage/"
                )
            )
            
        files.df <- df %>% 
            dplyr::group_split(bp.type)

        cor.mat <- lapply(files.df, function(f){
            res <- lapply(f$file.name, function(x){
                temp <- switch(action,
                    "zscore" = {
                        fread(x, select = "z")
                    },
                    "ratio" = {
                        dt <- fread(x, select = c("case", "control"))
                        norm.case <- dt$case/sum(dt$case, na.rm = TRUE)
                        norm.control <- dt$control/sum(dt$control, na.rm = TRUE)
                        ratio <- norm.case/norm.control
                    }
                )
                return(temp)
            })
            res <- do.call(cbind, res)
            colnames(res) <- paste0("V", 1:ncol(res))

            # remove any values of inf and NAs
            rows.to.keep <- which(is.finite(rowSums(res)))
            res <- res[rows.to.keep, ]

            # calculate correlation coefficients
            cor.coef <- cor(
                res, 
                use = "everything", 
                method = "pearson"
            )
            cor.coef <- as.numeric(cor.coef[1,])

            # keep non 1 correlation coefficients
            cor.coef <- cor.coef[cor.coef != 1]    

            return(list(cor.coef, rep(unique(f$bp.type), length(cor.coef))))
        })

        cors <- unlist(sapply(cor.mat, `[[`, 1))
        exps <- unlist(sapply(cor.mat, `[[`, 2))
        res <- data.table(exp = exps, cors = cors, kmer = kmer, action = action)

        return(res)
    })
    res <- rbindlist(res)
    return(res)
})
cor.mat.original <- rbindlist(cor.mat)

cor.mat <- as_tibble(cor.mat.original) %>% 
    dplyr::mutate(
        action = ifelse(
            action == "ratio", 
            "Prob. ratios", "Z-scores"
        ),
        action = factor(action, levels = c("Z-scores", "Prob. ratios")),
        kmer = paste0(kmer, "-mer"),
        kmer = as.factor(kmer)
    )
    
# combine both plots side-by-side
p_exactbreaks <- as_tibble(df_exactbreaks) %>%
    dplyr::mutate(
        overlaps = overlaps*100,
        group = "Exact breakage sites"
    ) %>%
    ggplot(aes(x = overlaps)) + 
    geom_density(fill = "#69b3a2", color = "black", alpha = 0.75) + 
    facet_wrap(vars(group), nrow = 1) + 
    scale_x_continuous(limits = c(0, 100)) +
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) +
    coord_cartesian(xlim = c(0, 100)) +
    labs(
        x = "Overlaps, %",
        y = "Density"
    )

p_kmerfragility <- cor.mat %>% 
    ggplot(aes(x = cors, fill = kmer)) + 
    geom_density(alpha = 0.5) + 
    facet_wrap(vars(action), nrow = 1) + 
    scale_x_continuous(limits = c(0, NA)) +
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) +
    coord_cartesian(xlim = c(0, 1)) + 
    labs(
        x = "Pearson correlation coefficient",
        y = ""
    )

pdf(
    file = "../figures/exact_breaks/BreakOverlaps_and_KmerFragilityCorr.pdf",
    height = 5, width = 12
)
gridExtra::grid.arrange(p_exactbreaks, p_kmerfragility, nrow = 1, widths = c(1, 2))
save.plot <- dev.off()