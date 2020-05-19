#' Metropolis-Hastings random walk proposal and update
#'
#' Sampler uses a mean-zero Gaussian random walk as a proposal distribution.
#' 
#' @param x0 Initial location of sampler
#' @param lp Function to evaluate log-posterior at model parameters
#' @param cov.chol Cholesky decomposition of proposal covariance
#'
#' @importFrom stats rnorm runif
#' 
# @example examples/mhrw.R
#'
mhrw = function(x0, lp, cov.chol) {

  # dimension of parameter vector
  n = length(x0)
  
  # random variates for proposal
  z = rnorm(n = n)
  
  # proposal
  x = x0 + cov.chol %*% z
  
  # metropolis ratio
  #   Recall: proposal is symmetric, so proposal density is not needed
  lR = lp(x) - lp(x0)
  
  # accept/reject
  accept = log(runif(1)) <= lR
  if(accept) {
    res = x
  } else {
    res = x0
  }
  
  # package results
  list(x = res, accepted = accept)
}