#' Use bridged sampling to impute a complete dive trajectory consistent with observations
#'
#' The sampling method is designed to sample many trajectories simultaneously, 
#' so has an extra level of approximation in the proposal distributions.
#' 
#' @param B single-step transition matrix, after uniformization
#' @param L likelihood matrix where each column is the probability distribution
#'  for the state at each of the discrete transitions
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