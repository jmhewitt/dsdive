#' Metropolis-Hastings update for subsets of dive model parameters
#' 
#' The Metropolis-Hastings update in \code{gibbs.mhga.dsdive} uses the
#' Gaussian approximation to the posterior mode, conditional on 1) the family of 
#' complete or imputed dive trajectories, and 2) a search direction in which to 
#' move, as the proposal distribution.
#'
#' 
gibbs.mhgdir.dsdive = function(x0, ld0, lp0, cfg, priors, verbose = FALSE) {
  
  if(is.infinite(ld0 + lp0)) {
    stop('Attempting to take MH step from a location with 0 mass.')
  }
  
  # generate proposal direction
  d = rnorm(n = length(x0))
  d = d / sqrt(sum(d^2))

  # determine posterior mode along slice
  o = optim(par = 0, fn = function(lambda) {
    # munge parameters
    params = x0 + lambda * d
    params.prop.list = params.toList(par = params, spec = priors)
    # compute log-density
    prop.ld = dsdive_ld(cfg = cfg, params = params.prop.list) +
      params.prop.list$logJ
    # compute log posterior
    r = prop.ld + 
      dsdive.prior(par = params.prop.list, spec = priors, log = TRUE)
    # diagnostics
    if(FALSE) {
      print(lambda)
      message(r)
    }
    # return
    r
  }, method = 'BFGS', hessian = TRUE, control = list(maxit = 1e3, fnscale = -1))
  
  if(is.na(o$hessian)) {
    prop.sd = .01
  } else {
    # extract proposal sd
    prop.sd = sqrt(-1/o$hessian)
  }
  
  # compute original lambda
  lambda0 = -o$par
  
  # propose new lambda and parameters
  lambda = -rnorm(n = 1, sd = prop.sd)
  x = x0 + (lambda0 + lambda) * d
  
  # backtransform parameters and compute proposed likelihood/prior
  params = params.toList(par = x, spec = priors)
  ld = dsdive_ld(cfg = cfg, params = params) + params$logJ
  lp = dsdive.prior(par = params, spec = priors, log = TRUE)
  
  # accept/reject
  lR = (ld + lp) - (ld0 + lp0) + dnorm(x = lambda0, sd = prop.sd, log = TRUE) - 
    dnorm(x = lambda, sd = prop.sd, log = TRUE)
  if(log(runif(1)) <= lR) {
    accepted = 'ACCEPTED'
  } else {
    accepted = ''
    x = x0
    ld = ld0
    lp = lp0
    params = params.toList(par = x0, spec = priors)
  }
  
  if(verbose) {
    message(paste('   lR:', round(lR, 2), accepted, sep = ' '))
  }
  
  # return new parameter and log-posterior components
  list(x = x, ld = ld, lp = lp, params = params)
}