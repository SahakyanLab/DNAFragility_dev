FitGMM <- function(dat, ind, sigma3=2){
  # Fit to a model consisting of a pair of Gaussians. Note that nls is kind
  # of fussy about choosing "reasonable" starting guesses for the parameters.
  # It's even fussier if you use the default algorithm (i.e., Gauss-Newton)
  x <- pull(dat, x)
  y <- pull(dat, y)

  if(ind == "kmer_6"){
    fit <- nls(
      y ~ (C1 * exp(-(x-mean1)**2/(2 * sigma1**2)) +
           C2 * exp(-(x-mean2)**2/(2 * sigma2**2)) + min(y)),
      data=dat,
      start=list(C1=10, mean1=0, sigma1=sd(x),
                 C2=10, mean2=0, sigma2=1),
      algorithm="port"
    )
  } else {
    try.nls <- function(sigma3){
      return(
        nls(
          y ~ (C1 * exp(-(x-mean1)**2/(2 * sigma1**2)) +
              C2 * exp(-(x-mean2)**2/(2 * sigma2**2)) +
              C3 * exp(-(x-mean3)**2/(2 * sigma3**2)) + min(y)),
          data=dat,
          start=list(C1=10, mean1=0, sigma1=sd(x),
                    C2=10, mean2=0, sigma2=1,
                    C3=10, mean3=0, sigma3=sigma3),
          algorithm="port"
        )
      )
    }

    sigma.value <- 2
    operation <- function(sigma3){
      withRestarts(
        tryCatch({
          try.nls(sigma3 = sigma3)
        },
        error = function(e){
          invokeRestart("retry")
        }),
        retry = function(){
          sigma.value <<- sigma3+20
          stopifnot(sigma3 > 0)
          operation(sigma3+20)
        }
      )
    }  

    operation(sigma3 = sigma.value)
    fit <- try.nls(sigma3 = sigma.value)
  }

  # extract each parameter
  summary.fit <- summary(fit)
  summary.fit.params <- summary.fit$parameters
  
  # fit each Gaussian curve separately to data
  draw.from.gaussian <- function(xs, C, MEAN, SD) {
    return(C * exp(-(xs-MEAN)**2/(2 * SD**2)))
  }
  
  # fit curve 1
  fit_1 = draw.from.gaussian(
    xs = x,
    C = summary.fit.params["C1", 1],
    MEAN = summary.fit.params["mean1", 1],
    SD = summary.fit.params["sigma1", 1]
  )
  fit_1 = fit_1 + min(y)

  # fit curve 2
  fit_2 = draw.from.gaussian(
    xs = x,
    C = summary.fit.params["C2", 1],
    MEAN = summary.fit.params["mean2", 1],
    SD = summary.fit.params["sigma2", 1]
  )
  fit_2 = fit_2 + min(y)

  if(ind == "kmer_6"){
    return(list(fit_1, fit_2))
  } else {
    # fit curve 3
    fit_3 = draw.from.gaussian(
      xs = x,
      C = summary.fit.params["C3", 1],
      MEAN = summary.fit.params["mean3", 1],
      SD = summary.fit.params["sigma3", 1]
    )
    fit_3 = fit_3 + min(y)

    return(list(fit_1, fit_2, fit_3))
  }
}