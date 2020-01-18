#' Use bridged sampling to impute a complete dive trajectory consistent with observations
#'
#' The sampling method is designed to sample many trajectories simultaneously, 
#' so has an extra level of approximation in the proposal distributions.
#' 
#' @param B single-step transition matrix, after uniformization
#' @param x0 state at which chain begins
#' @param xN state at which chain ends
#' @param N number of steps to take
#' 
#' @importFrom Matrix sparseVector
#' 
#' @example examples/ffbs.segment.R
#' 
#' @export
#'
ffbs.segment = function(B, x0, xN, N) {
  
  if(N==0) {
    
    if(x0==xN) {
      s = x0
    } else {
      stop('Bridging is not possible; N==0, and x0!=xN.')
    }
    
  } else if(N==1) {
    
    if(B[x0,xN] > 0) {
      s = c(x0, xN)
    } else {
      stop('Bridging is not possible; N==1, and B indicates P(xN | x0) = 0.')
    }
    
  } else {
    
    # state space size
    m = nrow(B)
    
    # initialize forward-filtering vectors
    a = vector('list', N+1)
    
    # encode initial state distribution
    a[[1]] = sparseVector(1, x0, m)
    
    # forward filter
    Bt = t(B)
    for(i in 2:length(a)) {
      a[[i]] = Bt %*% a[[i-1]]
      # standardize for numerical stability
      a[[i]] = a[[i]] / sum(a[[i]])
    }
    
    #
    # backward sample
    #
    
    s = numeric(N+1)
    s[1] = x0
    s[N+1] = xN
    
    for(t in N:2) {
      p = as.numeric(B[,s[t+1]] * a[[t]])
      s[t] = sample(x = 1:m, size = 1, prob = p)
    }
    
  }
  
  s
}