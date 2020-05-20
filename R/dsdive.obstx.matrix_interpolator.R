#' Interpolators for probability transition matrices of partially observed CTMC
#' 
#' Given a range of model parameters, this function will compute the 
#' probability transition matrices for all combinations of a Continuous time 
#' Markov chain (CTMC) that is observed once, and  then again at a time of 
#' \code{tstep} units of time later.  Then, thin plate spline interpolating 
#' functions will be built to allow probability transition matrix entries to 
#' be approximated at different values of the input parameters.
#' 
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param beta.seq the range of depth bin transition model parameters
#' @param lambda.seq the range of depth bin transition rate model parameters
#' @param s0 the stage for which to compute the transition matrix
#' @param tstep.seq Range of times between observations of the CTMC
#' @param m Smoothing polynomial degree; this is the \code{m} argument in the 
#'   smoothing function, \code{fields::Tps}.
#' @param verbose \code{TRUE} to output progress of code
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