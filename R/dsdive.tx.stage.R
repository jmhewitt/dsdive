#' Compute stage transition parameters for dive trajectories across discrete depths
#'
#' Computes the probability that a dive will transition between any of the 3 
#' stages, PRIMARY DESCENT (PD), INTERMEDIATE BEHAVIORS (IB), and 
#' PRIMARY ASCENT (PA).
#'   
#' @param t0 vector of times at which transition probs. should be computed
#' @param d0 vector of depth bins at which transition probs. should be computed
#' @param sub.tx length 2 vector that specifies the first depth bin at which 
#'   transitions to the IB stage can occur and the probability that such 
#'   a transition occurs at the next depth transition
#' @param surf.tx parameter that specifies the probability the trajectory will 
#'   transition to the PA stage at the next depth transition
#' @param t0.dive Time at which dive started
#' 
#' @example examples/dsdive.tx.stage.R
#' 
#' @export
#' 
dsdive.tx.stage = function(t0, d0, sub.tx, surf.tx, t0.dive) {
  
  # shift times so that they are centered with respect to dive start time
  t0 = t0 - t0.dive
  
  # initialize results
  res = matrix(data = 0, nrow = 3, ncol = length(t0))
  
  # compute probabilities
  
  # probability of transition to stage 2
  res[1, d0 >= sub.tx[1]] = sub.tx[2]
  
  # probability of transition to stage 3
  res[2,] = plogis(surf.tx[1] + surf.tx[2] * t0)
  
  # probability of transition out of stage 3 remains 0
   
  res
}