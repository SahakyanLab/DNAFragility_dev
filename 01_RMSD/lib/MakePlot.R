MakePlot <- function(dat, k, curve.vals, nr.of.curves){
  # extract standard deviation values to compute 95% confidence intervals
  st.devs <- unname(curve.vals[[length(curve.vals)-1]])

  CI.lst <- numeric(length = nr.of.curves*2)
  CI.lst[1] <- -1.96*st.devs[[1]]
  CI.lst[2] <- 1.96*st.devs[[1]]

  if(nr.of.curves > 1){
    CI.lst[3] <- -1.96*st.devs[[2]]
    CI.lst[4] <- 1.96*st.devs[[2]]
  } 

  if(nr.of.curves == 3){
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
          y = curvefits[[1]]
        ),
        stat = "identity",
        color = "red",
        alpha = 0.8,
        size = 1.2) + 
      geom_vline(xintercept = CI.lst[1], col = "red", linetype = "dashed") + 
      geom_vline(xintercept = CI.lst[2], col = "red", linetype = "dashed") +
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
        geom_vline(xintercept = CI.lst[3], col = "blue", linetype = "dashed") + 
        geom_vline(xintercept = CI.lst[4], col = "blue", linetype = "dashed") + 
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
            "Red CI:   ", formatC(signif(CI.lst[1], 3), 3), " / ", formatC(signif(CI.lst[2], 3), 3), "\n",
            "Blue CI:  ", formatC(signif(CI.lst[3], 3), 3), " / ", formatC(signif(CI.lst[4], 3), 3)
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
        geom_vline(xintercept = CI.lst[5], col = "darkgreen", linetype = "dashed") + 
        geom_vline(xintercept = CI.lst[6], col = "darkgreen", linetype = "dashed") +
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
              "Red CI:   ", formatC(signif(CI.lst[1], 3), 3), " / ", formatC(signif(CI.lst[2], 3), 3), "\n",
              "Green CI: ", formatC(signif(CI.lst[5], 3), 3), " / ", formatC(signif(CI.lst[6], 3), 3), "\n",
              "Blue CI:  ", formatC(signif(CI.lst[3], 3), 3), " / ", formatC(signif(CI.lst[4], 3), 3)
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
