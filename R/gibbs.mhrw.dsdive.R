#' Random walk Metropolis-Hastings update for subsets of dive model parameters
#'
#' 
gibbs.mhrw.dsdive = function(x0, ld0, lp0, sigma.chol, cfg, priors, 
                             ind = 1:length(x0), verbose = FALSE) {
  
  if(is.infinite(ld0 + lp0)) {
    stop('Attempting to take MH-RW step from a location with 0 mass.')
  }
  
  # initialize proposal
  x = x0
  
  # propose new parameters on unrestricted support
  x[ind] = x[ind] + sigma.chol %*% rnorm(n = nrow(sigma.chol))
  
  # backtransform parameters and compute proposed likelihood/prior
  params = params.toList(par = x, spec = priors)
  ld = dsdive_ld(cfg = cfg, params = params) + params$logJ
  lp = dsdive.prior(par = params, spec = priors, log = TRUE)
  
  # accept/reject
  lR = (ld + lp) - (ld0 + lp0)
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