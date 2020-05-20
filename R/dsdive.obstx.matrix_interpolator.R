#' Compute probability transition matrix for partially observed CTMC
#' 
#' Given model parameters, this function will compute the probability transition 
#' matrix for a Continuous time Markov chain (CTMC) that is observed once, and 
#' then again at a time of \code{tstep} units of time later.  If 
#' \code{include.raw==TRUE}, then raw components of the probability transition 
#' matrix will be returned so that the transition matrix can be computed for 
#' arbitrary timesteps.
#' 
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param beta the current depth bin transition model parameters
#' @param lambda the current depth bin transition rate model parameters
#' @param s0 the stage for which to compute the transition matrix
#' @param tstep Time between observations of the CTMC
#' @param include.raw \code{TRUE} to include raw components of the probability
#'   transition matrix so that the matrix can be computed for arbitrary 
#'   timesteps
#' @param delta If \code{delta>0}, then the observation matrix and raw 
#'   components computed will be for a transition matrix whose generator is 
#'   perturbed to allow much faster computation.
#'   
#' @example examples/dsdive.obstx.matrix_interpolator.R
#' 
#' @importFrom fields Tps
#' 
#' @export
#' 
dsdive.obstx.matrix_interpolator = function(depth.bins, beta.seq, lambda.seq, 
                                            s0, tstep.seq, m, verbose = FALSE) {
  
  # parameter grid
  param.grid = expand.grid(beta = beta.seq, lambda = lambda.seq, 
                         tstep = tstep.seq)
  
  if(verbose) {
    message(paste('Computing transition matrices at', nrow(param.grid), 
                  'parameter/timestep combinations.', sep = ' '))
  }
  
  # transition matrices for all parameter combinations
  P = apply(param.grid, 1, function(params) {
    dsdive.obstx.matrix(
      depth.bins = depth.bins, beta = c(params['beta'], .5), 
      lambda = c(params['lambda'], 1, 1), s0 = 1, tstep = params['tstep'], 
      include.raw = TRUE, delta = 0)
  })
  
  # number of entries in transition matrix
  n = length(P[[1]]$obstx.mat)
  
  if(verbose) {
    message('Building interpolators')
  }
  
  # smoothing fits to all matrix entries
  fits = lapply(1:n, function(i) {
    
    if(verbose) {
      message(paste('  ', i, sep = ''))
    }
    
    # all log probabilities 
    logp = sapply(P, function(pmat) log(pmat$obstx.mat[i]) )
    # replace 0-probabilities with an "equivalently" small log-value
    logp[is.infinite(logp)] = -1e3
    # spline fit
    suppressWarnings(Tps(x = param.grid, Y = logp, m = m, lambda = 0))
  })
  
  # simple index into smooths
  key = matrix(1:n, nrow = sqrt(n))
  
  # build interpolator
  res = function(beta, lambda, tstep, i, j, log = FALSE) {
    r = predict(fits[[key[i,j]]], data.frame(beta, lambda, tstep))
    if(log) {
      r
    } else {
      exp(r)
    }
  }
 
 res 
}