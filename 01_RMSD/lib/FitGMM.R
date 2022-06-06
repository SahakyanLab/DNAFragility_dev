FitGMM <- function(dat, ind, C.value=1, sigma3=2, nr.of.curves){
  # Fit to a model consisting of a pair of Gaussians. Note that nls is kind
  # of fussy about choosing "reasonable" starting guesses for the parameters.
  # It's even fussier if you use the default algorithm (i.e., Gauss-Newton)
  x <- pull(dat, x)
  y <- pull(dat, y)

  # initilise weights
  we <- rep(1, length(y))

  # weightings towards edge of curves
  we[1] <- max(y)*y[1]
  we[length(we)] <- max(y)*y[length(y)]

  # weightings towards origin of breakpoint
  # we[1:50] <- 2
  # we[(length(we)-50):length(we)] <- 2

  # error capture
  init.list <- vector(mode = "list", length = 4)

  # error handling
  # operation.cvalue <- function(C.value){
  #   withRestarts(
  #     tryCatch({
  #       try.nls(C.value = C.value)
  #     },
  #     error = function(e){
  #       invokeRestart("retry")
  #     }),
  #     retry = function(){
  #       Coef.value <<- C.value+11
  #       if(C.value > 1200) return(2)
  #       operation.cvalue(C.value+11)
  #     }
  #   )
  # }

  # operation.sigma <- function(sigma3){
  #   withRestarts(
  #     tryCatch({
  #       try.nls(sigma3 = sigma3)
  #     },
  #     error = function(e){
  #       invokeRestart("retry")
  #     }),
  #     retry = function(){
  #       sigma.value <<- sigma3+20
  #       if(sigma3 >= 1000) return(2)
  #       operation.sigma(sigma3+20)
  #     }
  #   )
  # }

  control1 <- nls.control(maxiter = 1000)
  if(nr.of.curves == 1){
    try.nls <- function(C.value, algo){
      if(algo){
        return(
          nls(
            y ~ (C1 * exp(-(x-0)**2/(2 * sigma1**2)) + min(y)),
            data=dat,
            start=list(C1=C.value, sigma1=sd(x)),
            lower=c(rep(0, 2)),
            upper=c(rep(10000, 2)),
            weights=we,
            algorithm="port",
            control=control1
          )
        )
      } else {
        return(
          nls(
            y ~ (C1 * exp(-(x-0)**2/(2 * sigma1**2)) + min(y)),
            data=dat,
            start=list(C1=C.value, sigma1=sd(x)),
            weights=we,
            control=control1
          )
        )
      }
    }

    fit <- tryCatch({
      try.nls(C.value = C.value, algo = TRUE)
    }, error = function(e) return(2) )

    init.list[[4]] <- fit

    if(class(fit) != "nls"){
      fit <- tryCatch({
        try.nls(C.value = C.value, algo = FALSE)
      }, error = function(e) return(2))

      init.list[[4]] <- fit

      if(class(fit) == "nls"){
        if(!all(summary(fit)$parameters[c("C1"), 1] > 0)){
          return(init.list)
        }
      } else {
        return(init.list)
      }
    }
  } else if(nr.of.curves == 2){
    try.nls <- function(sigma3, C.value, algo){
      if(algo){
        return(
          nls(
            y ~ (C1 * exp(-(x-0)**2/(2 * sigma1**2)) +
                 C2 * exp(-(x-0)**2/(2 * sigma2**2)) + min(y)),
            data=dat,
            start=list(C1=C.value, sigma1=sd(x),
                       C2=C.value, sigma2=sigma3),
            lower=c(rep(0, 4)),
            upper=c(rep(10000, 4)),
            weights=we,
            algorithm="port"
          )
        )
      } else {
        return(
          nls(
            y ~ (C1 * exp(-(x-0)**2/(2 * sigma1**2)) +
                 C2 * exp(-(x-0)**2/(2 * sigma2**2)) + min(y)),
            data=dat,
            start=list(C1=C.value, sigma1=sd(x),
                       C2=C.value, sigma2=sigma3),
            weights=we,
            control=control1
          )
        )
      }
    }

    fit <- tryCatch({
      try.nls(sigma3 = sigma3, C.value = C.value, algo = TRUE)
    }, error = function(e) return(2) )

    init.list[[4]] <- fit

    if(class(fit) != "nls"){
      fit <- tryCatch({
        try.nls(sigma3 = sigma3, C.value = C.value, algo = FALSE)
      }, error = function(e) return(2))

      init.list[[4]] <- fit

      if(class(fit) == "nls"){
        if(!all(summary(fit)$parameters[c("C1", "C2"), 1] > 0)){
          return(init.list)
        }
      } else {
        return(init.list)
      }
    }
  } else {
    try.nls <- function(sigma3, C.value, algo=TRUE){
      if(algo){
        return(
          nls(
            y ~ (C1 * exp(-(x-0)**2/(2 * sigma1**2)) +
                 C2 * exp(-(x-0)**2/(2 * sigma2**2)) +
                 C3 * exp(-(x-0)**2/(2 * sigma3**2)) + min(y)),
            data=dat,
            start=list(C1=C.value, sigma1=sd(x),
                       C2=C.value, sigma2=sigma3,
                       C3=C.value, sigma3=sigma3*2),
            lower=c(rep(0, 6)),
            upper=c(rep(10000, 6)),
            weights=we,
            algorithm="port"
          )
        )
      } else {
        return(
          nls(
            y ~ (C1 * exp(-(x-0)**2/(2 * sigma1**2)) +
                 C2 * exp(-(x-0)**2/(2 * sigma2**2)) +
                 C3 * exp(-(x-0)**2/(2 * sigma3**2)) + min(y)),
            data=dat,
            start=list(C1=C.value, sigma1=sd(x),
                       C2=C.value, sigma2=sigma3,
                       C3=C.value, sigma3=sigma3*2),
            weights=we,
            control=control1
          )
        )

      }
    }

    fit <- tryCatch({
      try.nls(sigma3 = sigma3, C.value = C.value, algo = TRUE)
    }, error = function(e) return(2) )

    init.list[[4]] <- fit

    if(class(fit) != "nls"){
      fit <- tryCatch({
        try.nls(sigma3 = sigma3, C.value = C.value, algo = FALSE)
      }, error = function(e) return(2))

      init.list[[4]] <- fit

      if(class(fit) == "nls"){
        if(!all(summary(fit)$parameters[c("C1", "C2", "C3"), 1] > 0)){
          return(init.list)
        }
      } else {
        return(init.list)
      }
    }
  }

  summary.fit <- tryCatch({
    summary(fit)
  }, error = function(e) return(2))

  if(class(summary.fit) != "summary.nls"){
    return(2)
  }

  # extract each parameter
  summary.fit.params <- summary.fit$parameters
  
  # fit each Gaussian curve separately to data
  draw.from.gaussian <- function(xs, C, SD, miny) {
    return(C * exp(-(xs-0)**2/(2 * SD**2)) + min(miny))
  }
  
  # fit curve 1
  fit_1 = draw.from.gaussian(
    xs = x,
    C = summary.fit.params["C1", 1],
    SD = summary.fit.params["sigma1", 1],
    miny = min(y)
  )
  fit_1 = fit_1

  if(nr.of.curves == 1){
    return(
      list(
        fit_1, 
        list("curve.one" = 1),
        init.list[length(init.list)][[1]]
      )
    )
  } else if(nr.of.curves > 1){
    # fit curve 2
    fit_2 = draw.from.gaussian(
      xs = x,
      C = summary.fit.params["C2", 1],
      SD = summary.fit.params["sigma2", 1],
      miny = min(y)
    )
    fit_2 = fit_2

    # integrate curve 1
    integral.curve.one <- tryCatch({
      integrate(
        curve.one, 
        lower = x[1], 
        upper = x[length(x)],
        C = summary.fit.params["C1", 1],
        SD = summary.fit.params["sigma1", 1]
      )$value
    }, error = function(e){
      integrate(
        curve.one, 
        lower = x[1], 
        upper = x[length(x)],
        C = summary.fit.params["C1", 1],
        SD = summary.fit.params["sigma1", 1],
        rel.tol = 1e-15
      )$value
    })

    # integrate curve 2
    integral.curve.two <- tryCatch({
      integrate(
        curve.two, 
        lower = x[1], 
        upper = x[length(x)],
        C = summary.fit.params["C2", 1],
        SD = summary.fit.params["sigma2", 1]
      )$value
    }, error = function(e){
      integrate(
        curve.two, 
        lower = x[1], 
        upper = x[length(x)],
        C = summary.fit.params["C2", 1],
        SD = summary.fit.params["sigma2", 1],
          rel.tol = 1e-15
      )$value
    })

    if(nr.of.curves == 2){
      # integrate linear combination of gaussian curves
      integral.all.curves <- tryCatch({
        integrate(
          all.curves, 
          lower = x[1], 
          upper = x[length(x)],
          nr.of.curves = 2,
          C1 = summary.fit.params["C1", 1],
          C2 = summary.fit.params["C2", 1],
          SD1 = summary.fit.params["sigma1", 1],
          SD2 = summary.fit.params["sigma2", 1]
        )$value
      },error = function(e){
        integrate(
          all.curves, 
          lower = x[1], 
          upper = x[length(x)],
          nr.of.curves = 2,
          C1 = summary.fit.params["C1", 1],
          C2 = summary.fit.params["C2", 1],
          SD1 = summary.fit.params["sigma1", 1],
          SD2 = summary.fit.params["sigma2", 1],
          rel.tol = 1e-15
        )$value
      })

      # percentage contribution of each curve
      curve.one.contribution <- integral.curve.one/integral.all.curves
      curve.two.contribution <- 1-curve.one.contribution

      return(
        list(
          fit_1, fit_2, 
          list("curve.one" = curve.one.contribution, 
               "curve.two" = curve.two.contribution),
          init.list[length(init.list)][[1]]
        )
      )
    } else {
      # fit curve 3
      fit_3 = draw.from.gaussian(
        xs = x,
        C = summary.fit.params["C3", 1],
        SD = summary.fit.params["sigma3", 1],
        miny = min(y)
      )
      fit_3 = fit_3

      # integrate curve 3
      integral.curve.three <- tryCatch({
        integrate(
          curve.three, 
          lower = x[1], 
          upper = x[length(x)],
          C = summary.fit.params["C3", 1],
          SD = summary.fit.params["sigma3", 1]
        )$value
      }, error = function(e){
        integrate(
          curve.three, 
          lower = x[1], 
          upper = x[length(x)],
          C = summary.fit.params["C3", 1],
          SD = summary.fit.params["sigma3", 1],
          rel.tol = 1e-15
        )$value
      })

      # integrate linear combination of gaussian curves
      integral.all.curves <- tryCatch({
        integrate(
          all.curves, 
          lower = x[1], 
          upper = x[length(x)],
          nr.of.curves = 3,
          C1 = summary.fit.params["C1", 1],
          C2 = summary.fit.params["C2", 1],
          C3 = summary.fit.params["C3", 1],
          SD1 = summary.fit.params["sigma1", 1],
          SD2 = summary.fit.params["sigma2", 1],
          SD3 = summary.fit.params["sigma3", 1]
        )$value
      }, error = function(e){
        integrate(
          all.curves, 
          lower = x[1], 
          upper = x[length(x)],
          nr.of.curves = 3,
          C1 = summary.fit.params["C1", 1],
          C2 = summary.fit.params["C2", 1],
          C3 = summary.fit.params["C3", 1],
          SD1 = summary.fit.params["sigma1", 1],
          SD2 = summary.fit.params["sigma2", 1],
          SD3 = summary.fit.params["sigma3", 1],
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
          init.list[length(init.list)][[1]]
        )
      )
    }
  }
}