#' Full conditional distribution for dive stages given data and parameters
#'
#' The 3 stage dive model has two breakpoints where the dive stage changes.  
#' Given data and model parameters, this function computes the conditional 
#' density of one of the breakpoints given the other.  The function is useful 
#' for constructing a Gibbs sampler that can sample over unknown stages of a 
#' dive when the dive depth bins and durations are known.  For example, this 
#' is helpful when using the \code{crawl}-based computational strategy to 
#' explore the posterior distribution for model parameters.
#' 
#' Furthermore, this full conditional distribution serves as the full 
#' conditional posterior distribution for stage breakpoints because the stage 
#' breakspoints are conditionally independent from the model parameters given 
#' the dive depth bins and durations.
#' 
#' @param breaks the two indices where a new dive stage is entered
#' @param fixed.ind  which of the stage break points in \code{breaks} should be 
#'   held fixed
#' @param beta \eqn{2 x 3} matrix in which each column contains the diving 
#'  preference and directional persistence parameters for the DIVING, SUBMERGED, 
#'  and SURFACING dive stages.
#' @param lambda length 3 vector that specifies the transition rate, 
#'   respectively in the DIVING, SUBMERGED, and SURFACING stages.
#' @param sub.tx length 2 vector that specifies the first depth bin at which 
#'   transitions to the SUBMERGED stage can occur and the probability that such 
#'   a transition occurs at the next depth transition
#' @param surf.tx parameter that specifies the probability the trajectory will 
#'   transition to the SURFACING stage at the next depth transition
#' @param depths record of depth bins the trajectory should visit
#' @param durations record of amount of time spent in each depth bin
#' @param times times at which the depth bins should be visited
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param t0.dive Time at which dive started
#' @param t.stage2 time at which second stage was entered
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#'   
#'
#' 
#' @example examples/gaussapprox.R
#' 
#'
gaussapprox = function(logf, init, method = 'BFGS', 
                       optim.output = FALSE, lower = -Inf, upper = Inf,
                       control = list(fnscale =  -1)) {
  
  # get parameter dimension
  p = length(init)
  
  # optimize surface
  o = optim(par = init, fn = logf, method = method, lower = lower, 
            upper = upper, control = control, hessian = TRUE)
  
  # cholesky for precision matrix
  prec.chol = chol(-o$hessian)
  
  # proposal function
  rgaussapprox = function(n) {
    # generate random variates
    z = matrix(rnorm(n = n*p), nrow = p)
    # add correlation structure
    y = backsolve(r = prec.chol, x = z)
    # add mean structure
    x = sweep(x = y, MARGIN = 1, STATS = o$par, FUN = '+')
    if(p==1) { as.numeric(x) } else { x }
  }
  
  # log-integration constant
  cst = - p * 0.9189385332 + sum(log(diag(prec.chol)))
  
  # density function
  dgaussapprox = function(x, log = FALSE) {
    y = prec.chol %*% (x - o$par)
    r = cst - .5 * apply(y, 2, function(y) sum(y^2))
    if(log) { r } else { exp(r) }
  }
   
  
  #
  # package results
  #
  
  res = list(rgaussapprox = rgaussapprox, dgaussapprox = dgaussapprox)
  
  if(optim.output) { 
    res$optim.output = o
  }
  
  res
}