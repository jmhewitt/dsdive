#' Forward-filtering backward-sampling algorithm for hidden markov models
#'
#' Implementation of forward-filtering backward-sampling algorithm described in 
#' Appendix A of Rao and Teh (2013).
#' 
#' @references Rao, Vinayak, and Yee Whye Teh. "Fast MCMC sampling for Markov 
#'   jump processes and extensions." The Journal of Machine Learning Research 
#'   14.1 (2013): 3295-3320.
#' 
#' @param B single-step transition matrix, after uniformization
#' @param L likelihood matrix where each column is the probability distribution
#'  for the state at each of the discrete transitions
#' @param a0 initial distribution of state vector
#' 
#' @importFrom Matrix sparseVector
#' 
#' @example examples/ffbs.R
#' 
#' @export
#'
ffbs = function(B, L, a0 = L[,1]) {
  # compute forward filtering vectors
  a = ff(B = B, L = L, a0 = a0)
  # execute backward sampling
  bs(a = a, B = B, L = L)
}