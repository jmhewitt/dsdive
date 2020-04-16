#' Compute probability transition matrix for embedded DTMCs
#' 
#' Computes the probability transition matrix for the Discrete time Markov 
#' chain (DTMC) embedded in a uniformized Continuous time Markov chain (CTMC).
#' See Rao and Teh (2013) for an applied definition and 
#' use of infinitesimal generator matrices.
#' 
#' @references Rao, Vinayak, and Yee Whye Teh. "Fast MCMC sampling for Markov 
#'   jump processes and extensions." The Journal of Machine Learning Research 
#'   14.1 (2013): 3295-3320.
#'   
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param beta Directional preference model parameters.  See 
#'   \code{dsdive.tx.params} for more details.
#' @param lambda Diving rate model parameters.  See 
#'   \code{dsdive.tx.params} for more details.
#' @param s0 dive stage for which matrix should be computed
#' @param rate.uniformized uniformization rate, for standardizing transition
#'   rates between states
#'   
#' @example examples/dsdive.tx.matrix.uniformized.R
#' 
#' @importFrom Matrix sparseMatrix
#' 
#' @export
#' 
dsdive.tx.matrix.uniformized = function(depth.bins, beta, lambda, s0,
                                        rate.uniformized) {
  
  n = nrow(depth.bins)

  # initialize storage for nonzero entries (overcommit space)
  nd = n * 3
  x = numeric(length = nd)
  im = numeric(length = nd)
  jm = numeric(length = nd)
  next.entry = 1

  # loop over depth bins
  for(i in 1:n) {
    
    # depth bin transition parameters
    p = dsdive.tx.params(depth.bins = depth.bins, d0 = i, s0 = s0, beta = beta, 
                         lambda = lambda)
    
    # self-transitions
    self.tx = ifelse(any(p$probs>0), 1 - p$rate / rate.uniformized, 1)
    if(self.tx > 0) {
      x[next.entry] = self.tx
      im[next.entry] = i
      jm[next.entry] = i
      next.entry = next.entry + 1
    }
    
    # transitions to new depths
    for(k in 1:length(p$labels)) {
      bin.ind = p$labels[k]
      prob = (1-self.tx) * p$probs[k]
      if(prob > 0) {
        x[next.entry] = prob
        im[next.entry] = i
        jm[next.entry] = bin.ind
        next.entry = next.entry + 1
      }
    }
  }
  
  # remove entries with null probabilities so that sparse matrix is compressed
  keep.inds = which(x > 0)
  x = x[keep.inds]
  im = im[keep.inds]
  jm = jm[keep.inds]
  
  # build and return matrix
  sparseMatrix(i = im, j = jm, x = x, dims = rep(n, 2))
}