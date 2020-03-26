#' Forward sampling algorithm for Hidden Markov Models
#'
#' Implementation of forward sampling algorithm described in Appendix A of 
#' Rao and Teh (2013).
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
#' @import Matrix
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
  N = ncol(L) - 1
  
  if(length(B) != N) {
    msg = paste('Number of transition matrices (', length(B), ')', 
      ' does not match implied number of transitions (', N, ')', sep ='')
    stop(msg)
  }
  
  # state space size
  m = ncol(L)
  
  # initialize forward-filtering vectors
  a = vector('list', N+1)
  
  # extract initial state distribution
  a[[1]] = a0
  
  # forward filter
  if(N>0) {
    for(i in 2:length(a)) {
      a[[i]] = t(B[[i-1]]) %*% (a[[i-1]] * L[,i-1])
      # standardize for numerical stability
      a[[i]] = a[[i]] / sum(a[[i]])
    }
  }
  
  # return list of forward-filtering distributions
  a
}