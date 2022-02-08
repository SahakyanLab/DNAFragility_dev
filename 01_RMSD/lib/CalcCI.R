CalcCI <- function(dat, xs){
  # find borders of each distribution based on 95% confidence interval
  MEAN = mean(xs)
  SD = sd(xs)
  error <- (qnorm(0.975)*SD)/sqrt(nrow(dat))
  left <- MEAN-error
  right <- MEAN+error
  
  # obtain x-value range for the cut-off
  lower.bound <- min(dat$x)+which(xs > left)[1]
  upper.bound <- max(dat$x)-which(xs > right)[1]

  # return as data frame
  df <- data.frame(
    lower.bound = lower.bound,
    upper.bound = upper.bound
  )

  return(df)
}