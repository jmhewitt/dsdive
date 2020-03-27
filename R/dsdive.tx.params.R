#' Compute transition parameters for a given depth bin
#'
#' Computes elements of a transition matrix when given model parameters and 
#' time/space locations for a dive model that has 3 stages, descent, 
#' bottom/foraging, and ascent.
#'   
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param d0 index of depth bin for which to compute transition parameters
#' @param s0 stage for which to compute parameters
#' @param beta vector with directional preference parameters for descent and  
#'   ascent stages.
#' @param lambda vector with dive rate parameters for the descent, bottom, and 
#'   ascent stages.
#' 
#' @example examples/dsdive.tx.params.R
#' 
#' @importFrom stats plogis
#' 
#' @export
#' 
dsdive.tx.params = function(depth.bins, d0, s0, beta, lambda) {
  
  num.depths = nrow(depth.bins)
  
  # extract transition rate
  rate = lambda[s0] / (2 * depth.bins[d0, 2])
  
  # define transition probabilities between dive bins
  if(d0 == 1) {
    
    # surface can only transition downward
    nbrs = 2
    
    # can only leave surface depth bin in stages 1 or 2
    probs = ifelse(s0 < 3, 1, 0)
    
  } else if(d0 == 2) {
    
    nbrs = c(1,3)
    
    # can only return to surface bin in stage 3
    probs.down = ifelse(s0==3, beta[2], 1)
    probs = c(1-probs.down, probs.down)
    
  } else if(d0 == num.depths) {
    
    # bottom bin is a boundary; transition is deterministic
    nbrs = num.depths - 1
    probs = 1
    
  } else {
    
    nbrs = d0 + c(-1, 1)
    probs.down = ifelse(s0==1, beta[1], ifelse(s0==2, .5, beta[2]))
    probs = c(1-probs.down, probs.down)
    
  }
  
  # package results
  list(rate = rate, probs = probs, labels = nbrs)
}