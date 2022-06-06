MakePlot <- function(dat, k, curve.vals, nr.of.curves){
  fit.1.lower.bound <- CalcCI(dat, curve.vals[[1]])$lower.bound
  fit.1.upper.bound <- CalcCI(dat, curve.vals[[1]])$upper.bound
  
  CI.lst <- numeric(length = nr.of.curves*2)
  CI.lst[1] <- fit.1.lower.bound
  CI.lst[2] <- fit.1.upper.bound

  if(nr.of.curves > 1){
    fit.2.lower.bound <- CalcCI(dat, curve.vals[[2]])$lower.bound
    fit.2.upper.bound <- CalcCI(dat, curve.vals[[2]])$upper.bound

    CI.lst[3] <- fit.2.lower.bound
    CI.lst[4] <- fit.2.upper.bound
  } 

  if(nr.of.curves == 3){
    fit.3.lower.bound <- CalcCI(dat, curve.vals[[3]])$lower.bound
    fit.3.upper.bound <- CalcCI(dat, curve.vals[[3]])$upper.bound

    # if too many curves fitted, confidence interval will be on edges
    check.if.too.many.curves.fitted <- c(
      fit.1.lower.bound, fit.1.upper.bound,
      fit.2.lower.bound, fit.2.upper.bound,
      fit.3.lower.bound, fit.3.upper.bound
    )

    # if(any(check.if.too.many.curves.fitted > 270 | check.if.too.many.curves.fitted < -270)){
    #   return(2)
    # }

    CI.lst[5] <- fit.3.lower.bound
    CI.lst[6] <- fit.3.upper.bound
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
      geom_vline(xintercept = fit.1.lower.bound, col = "red", linetype = "dashed") + 
      geom_vline(xintercept = fit.1.upper.bound, col = "red", linetype = "dashed") +
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
        geom_vline(xintercept = fit.2.lower.bound, col = "blue", linetype = "dashed") + 
        geom_vline(xintercept = fit.2.upper.bound, col = "blue", linetype = "dashed") + 
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
            "Red CI:   ", formatC(signif(fit.1.lower.bound, 3), 3), " / ", formatC(signif(fit.1.upper.bound, 3), 3), "\n",
            "Blue CI:  ", formatC(signif(fit.2.lower.bound, 3), 3), " / ", formatC(signif(fit.2.upper.bound, 3), 3)
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
        geom_vline(xintercept = fit.3.lower.bound, col = "darkgreen", linetype = "dashed") + 
        geom_vline(xintercept = fit.3.upper.bound, col = "darkgreen", linetype = "dashed") +
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
              "Red CI:   ", formatC(signif(fit.1.lower.bound, 3), 3), " / ", formatC(signif(fit.1.upper.bound, 3), 3), "\n",
              "Green CI: ", formatC(signif(fit.3.lower.bound, 3), 3), " / ", formatC(signif(fit.3.upper.bound, 3), 3), "\n",
              "Blue CI:  ", formatC(signif(fit.2.lower.bound, 3), 3), " / ", formatC(signif(fit.2.upper.bound, 3), 3)
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