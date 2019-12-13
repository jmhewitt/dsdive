#' Metropolis-Hastings update for subsets of dive model parameters
#' 
#' The Metropolis-Hastings update in \code{gibbs.mhga.dsdive} uses the
#' Gaussian approximation to the posterior mode, conditional on the family of 
#' complete or imputed dive trajectories, as the proposal distribution.
#'
#' @importFrom MHadaptive makePositiveDefinite
#' 
gibbs.mhga.dsdive = function(x0, ld0, lp0, cfg, priors, mu.vec, 
                             sigma.chol, verbose = FALSE) {
  
  if(is.infinite(ld0 + lp0)) {
    stop('Attempting to take MH step from a location with 0 mass.')
  }
  
  # generate proposal
  z = rnorm(n = nrow(sigma.chol))
  x = mu.vec + sigma.chol %*% z
  
  # backtransform parameters and compute proposed likelihood/prior
  params = params.toList(par = x, spec = priors)
  ld = dsdive_ld(cfg = cfg, params = params) + params$logJ
  lp = dsdive.prior(par = params, spec = priors, log = TRUE)
  
  z.new = forwardsolve(sigma.chol, x - mu.vec)
  z.old = forwardsolve(sigma.chol, x0 - mu.vec)
  
  # accept/reject
  lR = (ld + lp) - (ld0 + lp0) + 
    .5 * t(z.old) %*% z.old - .5 * t(z.new) %*% z.new
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