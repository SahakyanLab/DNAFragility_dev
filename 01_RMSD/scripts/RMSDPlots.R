# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.experiment <- as.character(args[2])
chromosome <- as.character(args[3])
category <- as.character(args[4])
auto.fit <- as.logical(args[5])

# my.path="/Volumes/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/01_RMSD/scripts/"
# # breakpoint.experiment="22-sBLISS/Colibactin_Ecoli_induced_DSBs/Caco-2_etoposide_rep1"
# breakpoint.experiment="10-DSBCapture/NHEK_DSBs"
# chromosome=1
# category="Biological"
# auto.fit=FALSE
setwd(my.path)

suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))

source("../lib/FitGMM.R")
source("../lib/IntegralFunctions.R")
source("../lib/CalcCI.R")
source("../lib/MakePlot.R")

# load data sets
files <- list.files(
  path = paste0("../data/", breakpoint.experiment),
  pattern = "^rmsd_kmer_"
)
file.name <- paste0("kmer_", seq(from = 2, to = 10, by = 2))
files <- str_sort(files, numeric = TRUE)

data.sets <- lapply(1:length(files), function(x){
  out <- readRDS(
    file = paste0("../data/", breakpoint.experiment,"/", files[x])
  )

  limits <- length(out)/2-1
  out <- as_tibble(out) %>% 
    mutate(
      kmer = as.factor(file.name[x]), 
      x = -limits:(length(out)-limits-1)
    ) %>% 
    dplyr::rename(y = value)
  
  return(out)
})

data.sets <- do.call(rbind, data.sets)

if(auto.fit){
  # Overall RMSD plots
  plots <- data.sets %>%
    ggplot(aes(x = x, y = y)) + 
    geom_line(size = 0.8) + 
    facet_wrap(~kmer, ncol = 2, scales = "free_y") + 
    theme_bw() + 
    labs(
      x = "Position away from breakpoint",
      y = "RMSD"
    )

  ggsave(
    filename = paste0("../figures/", breakpoint.experiment, 
                      "/chr", chromosome, "_sidebyside_RMSD.pdf"),
    plot = plots, 
    height = 7, width = 7
  )
}

# GMM plots
data.sets <- data.sets %>%
  mutate(across(where(is.factor), as.character)) %>% 
  filter(kmer != "kmer_2" & kmer != "kmer_10")

p <- vector(mode = "list", length = 3)

# if auto.fit, it tries to fit 3 curves by default
# if fails, try 2 curves. If fails, try 1 curve
# if no curves can be fit, it continues to next iteration
if(auto.fit){
  fits <- c("kmer_4" = 3, "kmer_6" = 3, "kmer_8" = 3)
} else {
  # extract optimal number of curves for each exp
  curve.nr <- fread("../../Raw_data/org_file.csv")
  fits <- as_tibble(curve.nr) %>% 
    mutate(
      exp = paste0(`Fragmentation type`, "/", `Experiment folder`), 
      .before = 1
    ) %>% 
    dplyr::filter(exp == breakpoint.experiment) %>% 
    dplyr::select(kmer_4, kmer_6, kmer_8) %>% 
    tidyr::gather(key = "kmer", value = "curves") %>% 
    pull(curves, kmer)
}

pb <- txtProgressBar(min = 1, max = length(fits), style = 3)
for(k in 1:length(fits)){
  setTxtProgressBar(pb, k)
  if(!is.na(fits[k])){
    dat <- data.sets %>% filter(kmer == names(fits[k]))
    curvefits <- FitGMM(dat = dat, ind = names(fits[k]), nr.of.curves = fits[k])

    if(auto.fit){
      if(class(curvefits[length(curvefits)][[1]]) != "nls"){
        fits[k] <- 2
        curvefits <- FitGMM(dat = dat, ind = names(fits[k]), nr.of.curves = fits[k])

        if(class(curvefits[length(curvefits)][[1]]) != "nls"){
          fits[k] <- 1
          curvefits <- FitGMM(dat = dat, ind = names(fits[k]), nr.of.curves = fits[k])
        }

        if(all(is.na(unlist(curvefits[length(curvefits)-2])))){
          p <- p[1:(length(p)-1)]
          break
        }
      }
    }

    output <- MakePlot(dat = dat, k = names(fits[k]), 
                       curve.vals = curvefits, 
                       nr.of.curves = fits[k])
    output.plot <- output[[1]]

    if(auto.fit){
      if(class(output.plot)[1] == "numeric"){
        fits[k] <- 2
        curvefits <- FitGMM(dat = dat, ind = names(fits[k]), nr.of.curves = fits[k])

        if(class(curvefits[length(curvefits)][[1]]) != "nls" & names(fits[k]) != "kmer_8"){
          fits[k] <- 1
          curvefits <- FitGMM(dat = dat, ind = names(fits[k]), nr.of.curves = fits[k])
        }

        if(all(is.na(unlist(curvefits[length(curvefits)-2])))){
          p <- p[1:(length(p)-1)]
          break
        }

        output <- MakePlot(dat = dat, k = names(fits[k]), 
                          curve.vals = curvefits, 
                          nr.of.curves = fits[k])
        output.plot <- output[[1]]
      } else {
        p[[k]] <- output.plot
      }
    } else {
      p[[k]] <- output.plot
    }

    # percent contribution of each gaussian curve towards breakage
    df <- as_tibble(curvefits[length(curvefits)-1][[1]])
    df[1,] <- df[1,]*100

    # range of influence based on 95 percent confidence intervals
    output.CIlst <- output[[2]]
    curves <- output.CIlst[seq(2, length(output.CIlst), 2)]-
              output.CIlst[seq(1, length(output.CIlst), 2)]

    # sd of curves
    SD <- sapply(1:(length(curvefits)-2), function(x){
      vals <- curvefits[[x]]

      left <- t.test(vals, conf.level = 0.68)$conf.int[1]
      right <- t.test(vals, conf.level = 0.68)$conf.int[2]

      lower.bound <- which(vals > left)[1]
      upper.bound <- which(vals > right)
      upper.bound <- upper.bound[length(upper.bound)]

      return((upper.bound-lower.bound)/2)
    })

    SD <- c(SD, rep(NA_real_, 3-length(SD)))

    # peak intensity of each gaussian curve
    peak.intensity <- sapply(1:fits[k], function(x){
      max(curvefits[[x]])
    })

    peak.intensity <- c(peak.intensity, rep(NA_real_, 3-length(peak.intensity)))

    # combine results
    df <- rbind(df, curves, SD, peak.intensity)
    df <- df %>% 
      mutate(
        exp = breakpoint.experiment,
        category = category,
        rowid = c("contribution", "ranges", "SD", "peak.intensty"), 
        .before = 1,
      )

    if(!auto.fit){
      # extract cut-off values for 
      cutoff.df <- fread("../data/Ranges_cutoffs.csv")
      cutoff.df <- as_tibble(cutoff.df) %>% 
        pull(Cluster, names(fits[k]))

      temp <- df %>% 
        filter(rowid == "ranges") %>% 
        dplyr::select(contains("curve")) %>% 
        tidyr::gather(key = "ranges", value = "curves") %>% 
        pull(ranges, curves)

      # short range
      short.range <- which(as.integer(names(temp)) <= as.integer(names(cutoff.df[3])))
      if(any(short.range)){
        if(length(short.range) > 1){
          # if multiple curves fall within the same range
          # then at least 1 is a redundant curve to be removed
          # after this, update the curve number on org_file.csv
          contribution.temp <- df %>% 
            filter(rowid == "contribution") %>% 
            dplyr::select(-c(1:3)) %>% 
            tidyr::gather(key = "curves", value = "contribution") %>% 
            pull(contribution, curves)

          max.index <- which.max(contribution.temp[short.range])
          temp[short.range[max.index]] <- cutoff.df[short.range[max.index]]
          colnames(df)[4:6][short.range[max.index]] <- cutoff.df[short.range[max.index]]
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
            filter(rowid == "contribution") %>% 
            dplyr::select(-c(1:3)) %>% 
            tidyr::gather(key = "curves", value = "contribution") %>% 
            pull(contribution, curves)

          max.index <- which.max(contribution.temp[mid.range])
          temp[mid.range[max.index]] <- cutoff.df[mid.range[max.index]]
          colnames(df)[4:6][mid.range[max.index]] <- cutoff.df[mid.range[max.index]]
        } else {
          temp[mid.range] <- cutoff.df[2]
          colnames(df)[4:6][mid.range] <- cutoff.df[2]
        }
      }

      # long range
      long.range <- which(as.integer(names(temp)) > as.integer(names(cutoff.df[2])))
      if(any(long.range)){
        if(length(long.range) > 1){
          # if multiple curves fall within the same range
          # then at least 1 is a redundant curve to be removed
          # after this, update the curve number on org_file.csv
          contribution.temp <- df %>% 
            filter(rowid == "contribution") %>% 
            dplyr::select(-c(1:3)) %>% 
            tidyr::gather(key = "curves", value = "contribution") %>% 
            pull(contribution, curves)

          max.index <- which.max(contribution.temp[long.range])
          temp[long.range[max.index]] <- cutoff.df[long.range[max.index]]
          colnames(df)[4:6][long.range[max.index]] <- cutoff.df[long.range[max.index]]
        } else {
          temp[long.range] <- "long.range"
          colnames(df)[4:6][long.range] <- "long.range"
        }
      }

      # drop any non-range-assigned curves
      to.remove <- 3+which(c("curve.one", "curve.two", "curve.three") %in% unname(temp))
      if(length(to.remove) > 0){
        df <- df %>% 
          dplyr::select(-all_of(to.remove))

        # update the curve number on curve_counts and org_file csv files
        curve.counts.csv <- fread(
          file = paste0("../data/", names(fits[k]),"_curve_counts.csv")
        )
        row.to.replace <- which(curve.counts.csv$exp == breakpoint.experiment)
        updated.curve.numbers <- ncol(df[4:length(df)])
        curve.counts.csv[row.to.replace, "Count"] <- updated.curve.numbers
        fwrite(
          x = curve.counts.csv,
          file = paste0("../data/", names(fits[k]),"_curve_counts.csv")
        )

        curve.nr[, exp := paste0(`Fragmentation type`, "/", `Experiment folder`)]
        row.to.replace <- which(curve.nr$exp == breakpoint.experiment)
        curve.nr[row.to.replace, (names(fits[k])) := updated.curve.numbers]
        curve.nr[, exp := NULL]
        fwrite(
            x = curve.nr, 
            file = "../../Raw_data/org_file.csv"
        )
      }
    }

    fwrite(
      x = df,
      file = paste0("../data/", breakpoint.experiment, 
                    ifelse(auto.fit, "/key", "/new"),
                    "_stats_", names(fits[k]), ".csv"),
      row.names = FALSE
    )
  }
  # if(k == 1) print(p) # temp. avoid an unsolved ggplot bug
}
close(pb)
 
# repeat process to get updated curve statistics
if(!auto.fit){
  # extract optimal number of curves for each exp
  curve.nr <- fread("../../Raw_data/org_file.csv")
  fits <- as_tibble(curve.nr) %>% 
    mutate(
      exp = paste0(`Fragmentation type`, "/", `Experiment folder`), 
      .before = 1
    ) %>% 
    dplyr::filter(exp == breakpoint.experiment) %>% 
    dplyr::select(kmer_4, kmer_6, kmer_8) %>% 
    tidyr::gather(key = "kmer", value = "curves") %>% 
    pull(curves, kmer)

  pb <- txtProgressBar(min = 1, max = length(fits), style = 3)
  for(k in 1:length(fits)){
    setTxtProgressBar(pb, k)
    if(!is.na(fits[k])){
      dat <- data.sets %>% filter(kmer == names(fits[k]))
      curvefits <- FitGMM(dat = dat, ind = names(fits[k]), nr.of.curves = fits[k])
      output <- MakePlot(dat = dat, k = names(fits[k]), 
                         curve.vals = curvefits, 
                         nr.of.curves = fits[k])
      output.plot <- output[[1]]
      p[[k]] <- output.plot

      # percent contribution of each gaussian curve towards breakage
      df <- as_tibble(curvefits[length(curvefits)-1][[1]])
      df[1,] <- df[1,]*100

      # range of influence based on 95 percent confidence intervals
      output.CIlst <- output[[2]]
      curves <- output.CIlst[seq(2, length(output.CIlst), 2)]-
                output.CIlst[seq(1, length(output.CIlst), 2)]

      # sd of curves
      SD <- sapply(1:(length(curvefits)-2), function(x){
        vals <- curvefits[[x]]

        left <- t.test(vals, conf.level = 0.68)$conf.int[1]
        right <- t.test(vals, conf.level = 0.68)$conf.int[2]

        lower.bound <- which(vals > left)[1]
        upper.bound <- which(vals > right)
        upper.bound <- upper.bound[length(upper.bound)]

        return((upper.bound-lower.bound)/2)
      })

      SD <- c(SD, rep(NA_real_, 3-length(SD)))

      # peak intensity of each gaussian curve
      peak.intensity <- sapply(1:fits[k], function(x){
        max(curvefits[[x]])
      })

      peak.intensity <- c(peak.intensity, rep(NA_real_, 3-length(peak.intensity)))

      # combine results
      df <- rbind(df, curves, SD, peak.intensity)
      df <- df %>% 
        mutate(
          exp = breakpoint.experiment,
          category = category,
          rowid = c("contribution", "ranges", "SD", "peak.intensty"), 
          .before = 1,
        )

      if(!auto.fit){
        # extract cut-off values for 
        cutoff.df <- fread("../data/Ranges_cutoffs.csv")
        cutoff.df <- as_tibble(cutoff.df) %>% 
          pull(Cluster, names(fits[k]))

        temp <- df %>% 
          filter(rowid == "ranges") %>% 
          dplyr::select(contains("curve")) %>% 
          tidyr::gather(key = "ranges", value = "curves") %>% 
          pull(ranges, curves)

        # short range
        short.range <- which(as.integer(names(temp)) <= as.integer(names(cutoff.df[3])))
        if(any(short.range)){
          if(length(short.range) > 1){
            # if multiple curves fall within the same range
            # then at least 1 is a redundant curve to be removed
            # after this, update the curve number on org_file.csv
            contribution.temp <- df %>% 
              filter(rowid == "contribution") %>% 
              dplyr::select(-c(1:3)) %>% 
              tidyr::gather(key = "curves", value = "contribution") %>% 
              pull(contribution, curves)

            max.index <- which.max(contribution.temp[short.range])
            temp[short.range[max.index]] <- cutoff.df[short.range[max.index]]
            colnames(df)[4:6][short.range[max.index]] <- cutoff.df[short.range[max.index]]
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
              filter(rowid == "contribution") %>% 
              dplyr::select(-c(1:3)) %>% 
              tidyr::gather(key = "curves", value = "contribution") %>% 
              pull(contribution, curves)

            max.index <- which.max(contribution.temp[mid.range])
            temp[mid.range[max.index]] <- cutoff.df[mid.range[max.index]]
            colnames(df)[4:6][mid.range[max.index]] <- cutoff.df[mid.range[max.index]]
          } else {
            temp[mid.range] <- cutoff.df[2]
            colnames(df)[4:6][mid.range] <- cutoff.df[2]
          }
        }

        # long range
        long.range <- which(as.integer(names(temp)) > as.integer(names(cutoff.df[2])))
        if(any(long.range)){
          if(length(long.range) > 1){
            # if multiple curves fall within the same range
            # then at least 1 is a redundant curve to be removed
            # after this, update the curve number on org_file.csv
            contribution.temp <- df %>% 
              filter(rowid == "contribution") %>% 
              dplyr::select(-c(1:3)) %>% 
              tidyr::gather(key = "curves", value = "contribution") %>% 
              pull(contribution, curves)

            max.index <- which.max(contribution.temp[long.range])
            temp[long.range[max.index]] <- cutoff.df[long.range[max.index]]
            colnames(df)[4:6][long.range[max.index]] <- cutoff.df[long.range[max.index]]
          } else {
            temp[long.range] <- "long.range"
            colnames(df)[4:6][long.range] <- "long.range"
          }
        }

        # drop any non-range-assigned curves
        to.remove <- 3+which(c("curve.one", "curve.two", "curve.three") %in% unname(temp))
        if(length(to.remove) > 0){
          df <- df %>% 
            dplyr::select(-all_of(to.remove))

          # update the curve number on curve_counts and org_file csv files
          curve.counts.csv <- fread(
            file = paste0("../data/", names(fits[k]),"_curve_counts.csv")
          )
          row.to.replace <- which(curve.counts.csv$exp == breakpoint.experiment)
          updated.curve.numbers <- ncol(df[4:length(df)])
          curve.counts.csv[row.to.replace, "Count"] <- updated.curve.numbers
          fwrite(
            x = curve.counts.csv,
            file = paste0("../data/", names(fits[k]),"_curve_counts.csv")
          )

          curve.nr[, exp := paste0(`Fragmentation type`, "/", `Experiment folder`)]
          row.to.replace <- which(curve.nr$exp == breakpoint.experiment)
          curve.nr[row.to.replace, (names(fits[k])) := updated.curve.numbers]
          curve.nr[, exp := NULL]
          fwrite(
              x = curve.nr, 
              file = "../../Raw_data/org_file.csv"
          )

          # re-fit last curve again for updated stats and other metrics
          # extract optimal number of curves for each exp
          curve.nr <- fread("../../Raw_data/org_file.csv")
          fits <- as_tibble(curve.nr) %>% 
            mutate(
              exp = paste0(`Fragmentation type`, "/", `Experiment folder`), 
              .before = 1
            ) %>% 
            dplyr::filter(exp == breakpoint.experiment) %>% 
            dplyr::select(kmer_4, kmer_6, kmer_8) %>% 
            tidyr::gather(key = "kmer", value = "curves") %>% 
            pull(curves, kmer)

          dat <- data.sets %>% filter(kmer == names(fits[k]))
          curvefits <- FitGMM(dat = dat, ind = names(fits[k]), nr.of.curves = fits[k])
          output <- MakePlot(dat = dat, k = names(fits[k]), 
                            curve.vals = curvefits, 
                            nr.of.curves = fits[k])
          output.plot <- output[[1]]
          p[[k]] <- output.plot

          # percent contribution of each gaussian curve towards breakage
          df <- as_tibble(curvefits[length(curvefits)-1][[1]])
          df[1,] <- df[1,]*100

          # range of influence based on 95 percent confidence intervals
          output.CIlst <- output[[2]]
          curves <- output.CIlst[seq(2, length(output.CIlst), 2)]-
                    output.CIlst[seq(1, length(output.CIlst), 2)]

          # sd of curves
          SD <- sapply(1:(length(curvefits)-2), function(x){
            vals <- curvefits[[x]]

            left <- t.test(vals, conf.level = 0.68)$conf.int[1]
            right <- t.test(vals, conf.level = 0.68)$conf.int[2]

            lower.bound <- which(vals > left)[1]
            upper.bound <- which(vals > right)
            upper.bound <- upper.bound[length(upper.bound)]

            return((upper.bound-lower.bound)/2)
          })

          SD <- c(SD, rep(NA_real_, 3-length(SD)))

          # peak intensity of each gaussian curve
          peak.intensity <- sapply(1:fits[k], function(x){
            max(curvefits[[x]])
          })

          peak.intensity <- c(peak.intensity, rep(NA_real_, 3-length(peak.intensity)))

          # combine results
          df <- rbind(df, curves, SD, peak.intensity)
          df <- df %>% 
            mutate(
              exp = breakpoint.experiment,
              category = category,
              rowid = c("contribution", "ranges", "SD", "peak.intensty"), 
              .before = 1,
            )

          # extract cut-off values for 
          cutoff.df <- fread("../data/Ranges_cutoffs.csv")
          cutoff.df <- as_tibble(cutoff.df) %>% 
            pull(Cluster, names(fits[k]))

          temp <- df %>% 
            filter(rowid == "ranges") %>% 
            dplyr::select(contains("curve")) %>% 
            tidyr::gather(key = "ranges", value = "curves") %>% 
            pull(ranges, curves)

          # short range
          short.range <- which(as.integer(names(temp)) <= as.integer(names(cutoff.df[3])))
          if(any(short.range)){
            temp[short.range] <- cutoff.df[3]
            colnames(df)[4:6][short.range] <- cutoff.df[3]
          }

          # mid range
          mid.range <- which(
            (as.integer(names(temp)) > as.integer(names(cutoff.df[3]))) & 
            (as.integer(names(temp)) <= as.integer(names(cutoff.df[2])))
          )
          if(any(mid.range)){
            temp[mid.range] <- cutoff.df[2]
            colnames(df)[4:6][mid.range] <- cutoff.df[2]
          }

          # long range
          long.range <- which(as.integer(names(temp)) > as.integer(names(cutoff.df[2])))
          if(any(long.range)){
            temp[long.range] <- "long.range"
            colnames(df)[4:6][long.range] <- "long.range"
          }
        }
      }

      fwrite(
        x = df,
        file = paste0("../data/", breakpoint.experiment, 
                      ifelse(auto.fit, "/key", "/new"),
                      "_stats_", names(fits[k]), ".csv"),
        row.names = FALSE
      )
    }
    # if(k == 1) print(p) # temp. avoid an unsolved ggplot bug
  }
  close(pb)
}

ggsave(
  filename = paste0("../figures/", breakpoint.experiment, 
                    "/chr", chromosome, "_RMSD_all_",
                    ifelse(auto.fit, "GMMFit.pdf", "GMMFit_optimised.pdf")),
  plot = cowplot::plot_grid(plotlist = p, axis = "b", ncol = 1), 
  height = 15, width = 7
)