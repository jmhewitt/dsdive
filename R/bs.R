#' Backwards sampling algorithm for Hidden Markov Models
#' 
#' Implementation of backwards sampling algorithm described in Appendix A of 
#' Rao and Teh (2013).
#' 
#' @references Rao, Vinayak, and Yee Whye Teh. "Fast MCMC sampling for Markov 
#'   jump processes and extensions." The Journal of Machine Learning Research 
#'   14.1 (2013): 3295-3320.
#'   
#' @param a collection of forward filtering distributions
#' @param B collection of uniformized single-step transition matrices
#' @param L likelihood matrix where each column is the probability distribution
#'  for the state at each of the discrete transitions
#' 
#' @importFrom Matrix sparseVector
#' 
#' @example examples/ffbs.R
#' 
#' @export
#'
bs = function(a, B, L) {

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
  
  if(length(a) != N+1) {
    stop('a must include state distributions from inital position to end.')
  }
  
  # state space size
  m = nrow(L)
  
  #
  # backward sample
  #
  
  s = numeric(N+1)
  
  p = as.numeric(L[,N+1] * a[[N+1]])
  s[N+1] = sample(x = 1:m, size = 1, prob = p)
  
  if(N>0) {
    for(t in N:1) {
      p = as.numeric(B[[t]][,s[t+1]] * a[[t]] * L[,t])
      s[t] = sample(x = 1:m, size = 1, prob = p)
    }
  }
  
  s
}