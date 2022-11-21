CalcCI <- function(dat, xs){
  # find borders of each distribution based on 95% confidence interval
  MEAN = mean(xs)
  SD = sd(xs)

  # qnorm is the quantile function for the normal distribution. 
  # You pick 0.975 to get a two-sided confidence interval. 
  # This gives 2.5% of the probability in the upper tail and 2.5% in the lower tail.

  # The below "error" calculation is not itself a confidence interval, 
  # it computes the half-width of the interval.
  # You need two numbers (left and right end or mean and half width both work).
  error <- (qnorm(p=0.975)*SD)/sqrt(nrow(dat))
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