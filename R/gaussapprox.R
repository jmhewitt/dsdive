#' Compute a Gaussian approximation to a density
#'
#' Computes a Gaussian approximation to the target density \code{logf}, and 
#' additionally return methods to aid in sampling from the approximation.  The 
#' approximation is computed numerically via the \code{stats::optim} routine.
#' 
#' @param logf log of the density for which an approximation should be built.
#' @param init initial guess for mode of \code{logf}
#' @param method method to be used in \code{stats::optim}
#' @param lower lower bound value to pass to \code{stats::optim}
#' @param upper upper bound value to pass to \code{stats::optim}
#' @param control control list to pass to \code{stats::optim}
#' @param optim.output \code{TRUE} to return output from \code{stats::optim}
#' 
#' @importFrom stats optim rnorm 
#' 
# @example examples/gaussapprox.R
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
  
  class(res) = 'gaussapprox'
  
  if(optim.output) { 
    res$optim.output = o
  }
  
  res
}