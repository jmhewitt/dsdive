#' Basic implementation of stochastic optimization when the target function is only approximately known
#' 
#' Experimental method because its theoretical convergence and recommendations 
#' for setting parameters has not been fully explored.
#' 
#' @param par initial parameter values for function
#' @param fn a function to be maximized
#' @param maxit maximum number of function evaluations before failing to converge
#' @param a decay parameter for step size
#' @param c decay parameter for step size
#' @param A decay parameter for step size
#' @param alpha decay parameter for step size
#' @param gamma decay parameter for step size
#' @param tol Relative convergence criterion
#' 
#' @example examples/spsa.R
#' 
#' 
#'
spsa = function(par, fn, maxit = 1e2, a = 1, c = 1, A = 1, alpha = .5, 
                gamma = 1, tol = 1e-3) {
  
  # get model dimension
  p = length(par)
  
  # optimize
  err = tol
  converged = 1
  theta = par
  ghat.last = rep(1e-6, p)
  for(k in 1:maxit) {
    
    # evaluate gain sequence
    ak = a/(A+k)^alpha
    ck = c/k^gamma
    
    # random perturbation
    delta = 2 * round(runif(p)) - 1
    
    # approximate gradient, and update parameter if finite gradient estimate
    fplus = fn(theta + ck*delta)
    if(is.finite(fplus)) {
      
      fminus = fn(theta - ck*delta)
      if(is.finite(fminus)) {
        
        # update parameter
        ghat = (fplus - fminus) / (2 * ck * delta)
        theta = theta - ak * ghat
        
        # check convergence
        err = sqrt(sum((ghat - ghat.last)^2) / sum(ghat.last^2))
        ghat.last = ghat
        if(err < tol) {
          converged = 0
          break
        }
      }
    }
  }
  
  # package results
  list(
    theta = theta,
    converged = converged
  )
}