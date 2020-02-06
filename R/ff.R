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
ff = function(B, L, a0 = L[,1]) {

  if(!inherits(B, 'list')) {
    stop('B must be a list of single-step transition matrices')
  }
  
  # number of transitions to sample
  N = ncol(L)
  
  if(length(B) != N) {
    msg = paste('Number of transition matrices (', length(B), ')', 
      ' does not match implied number of transitions (', N, ')', sep ='')
    stop(msg)
  }
  
  # state space size
  m = nrow(B[[1]])
  
  # initialize forward-filtering vectors
  a = vector('list', N+1)
  
  # extract initial state distribution
  a[[1]] = a0
  
  # forward filter
  for(i in 2:length(a)) {
    a[[i]] = t(B[[i-1]]) %*% (a[[i-1]] * L[,i-1])
    # standardize for numerical stability
    a[[i]] = a[[i]] / sum(a[[i]])
  }
  
  # return list of forward-filtering distributions
  a
}