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
#' @param cl cluster over which computations can be parallelized
#'   
#' @example examples/dsdive.obstx.matrix_interpolator.R
#' 
#' @importFrom fields Tps
#' @importFrom parallel parLapply
#' 
#' @export
#' 
dsdive.obstx.matrix_interpolator = function(depth.bins, beta.seq, lambda.seq, 
                                            s0, tstep.seq, m, verbose = FALSE, 
                                            cl = NULL) {
  
  if(is.null(cl)) {
    lapply_impl = lapply
  } else {
    lapply_impl = function(X, FUN, ...) {
      parLapply(cl = cl, X = X, fun = FUN, ...)
    }
    clusterEvalQ(cl = cl, require(fields))
  }

  # parameter grid
  if(s0 != 2) {
    param.grid = expand.grid(beta = beta.seq, lambda = lambda.seq, 
                             tstep = tstep.seq)
  } else {
    param.grid = expand.grid(lambda = lambda.seq, tstep = tstep.seq)
  }
  
  
  if(verbose) {
    message(paste('Computing transition matrices at', nrow(param.grid), 
                  'parameter/timestep combinations.', sep = ' '))
  }
  
  # transition matrices for all parameter combinations
  P = lapply_impl(1:nrow(param.grid), function(i, param.grid, depth.bins, s0) {
    params = param.grid[i,]
    if(s0 == 1) {
      beta = c(params$beta, .5)
      lambda = c(params$lambda, 1, 1)
    } else if(s0 == 2) {
      beta = c(.5, .5)
      lambda = c(1, params$lambda, 1)
    } else if(s0 == 3) {
      beta = c(.5, params$beta)
      lambda = c(1, 1, params$lambda)
    }
    dsdive.obstx.matrix(
      depth.bins = depth.bins, beta = beta, lambda = lambda, s0 = s0, 
      tstep = params$tstep, include.raw = FALSE, delta = 0)
  }, param.grid, depth.bins, s0)
  
  # number of entries in transition matrix
  n = length(P[[1]])
  
  if(verbose) {
    message(paste('Building', n, 'interpolators.', sep = ' '))
  }
  
  # smoothing fits to all matrix entries
  fits = lapply_impl(1:n, function(i, P, param.grid, m) {
    
    if(verbose) {
      message(paste('  ', i, sep = ''))
    }
    
    # all log probabilities 
    logp = sapply(P, function(pmat) log(pmat[i]) )
    # replace 0-probabilities with an "equivalently" small log-value
    logp[is.infinite(logp)] = -1e3
    # spline fit
    tps.fit = suppressWarnings(Tps(x = param.grid, Y = logp, m = m, lambda = 0))
    # remove some large components not required for our predictions
    tps.fit$call = NULL
    tps.fit$matrices = NULL
    tps.fit$fitted.values = NULL
    tps.fit$fitted.values.null = NULL
    tps.fit$residuals = NULL
    tps.fit$gcv.grid = NULL
    tps.fit$warningTable = NULL
    tps.fit$uniquerows = NULL
    tps.fit$rep.info = NULL
    # return fit object
    tps.fit
  }, P, param.grid, m)
  
  # simple index into smooths
  key = matrix(1:n, nrow = sqrt(n))
  
  # build interpolator
  if(s0 != 2) {
    res = function(beta, lambda, tstep, i, j, log = FALSE) {
      r = predict(fits[[key[i,j]]], data.frame(beta, lambda, tstep))
      if(log) {
        r
      } else {
        exp(r)
      }
    }
  } else {
    res = function(beta, lambda, tstep, i, j, log = FALSE) {
      r = predict(fits[[key[i,j]]], data.frame(lambda, tstep))
      if(log) {
        r
      } else {
        exp(r)
      }
    }
  }
 
 res 
}