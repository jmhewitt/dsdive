#' Use bridged sampling to impute a complete dive trajectory consistent with observations
#'
#' The sampling method is designed to sample many trajectories simultaneously, 
#' so has an extra level of approximation in the proposal distributions.
#' 
#' @param B single-step transition matrix, after uniformization
#' @param x0 state at which chain begins
#' @param xN state at which chain ends
#' @param N maximum number of steps to consider
#' @param rate.uniformized uniformized poisson process transition rate
#' @param t time period between x0 and xN observations
#' 
#' @importFrom Matrix sparseVector
#' 
#' @example examples/dN.bridged.R
#' 
#' @export
#'
dN.bridged = function(B, x0, xN, N.max, rate.uniformized, t, log = FALSE) {
  
  # specify poisson distribution rate parameter
  lambdaT = rate.uniformized * t
  
  # build on top of poisson distribution component
  res = dpois(x = 0:N.max, lambda = lambdaT, log = TRUE)
  
  if(x0!=xN) {
    # impossible to get to endpoint without any transitions
    res[1] = -Inf
  }
  
  # initial state distribution
  a = sparseVector(1, x0, nrow(B))
  
  Bt = t(B)
  for(i in 2:length(res)) {
    
    # state dist'n. after one more step; standardize for numerical stability
    a = Bt %*% a
    a = a / sum(a)
    
    # add likelihood that the endpoint will be reached in i-1 transitions
    res[i] = res[i] + log(a[xN])
  }
  
  if(log) {
    res
  } else {
    res = exp(res)
    res / sum(res)
  }
}