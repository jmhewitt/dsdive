#' Conditional likelihood for completely observed dive trajectories
#' 
#' Computes the likelihood for depth bin durations and 
#' transitions in a completely observed dive trajectory conditional on stage 
#' transition times.  The conditioning treats the stage transition times as 
#' known quantities, rather than stochastic quantities.
#'
#' @param depths Depth bin indices in which the trajectory was observed
#' @param durations Amount of time spent in each depth bin
#' @param times Times at which each depth bin was entered
#' @param stages Record of dive stages at each depth bin
#' @param beta \eqn{2 x 3} matrix in which each column contains the diving 
#'  preference and directional persistence parameters for the DIVING, SUBMERGED, 
#'  and SURFACING dive stages.
#' @param lambda length 3 vector that specifies the transition rate, 
#'   respectively in the DIVING, SUBMERGED, and SURFACING stages.
#' @param sub.tx length 2 vector that specifies the first depth bin at which 
#'   transitions to the SUBMERGED stage can occur and the probability that such 
#'   a transition occurs at the next depth transition
#' @param surf.tx parameter that specifies the probability the trajectory will 
#'   transition to the SURFACING stage at the next depth transition
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param t0.dive Time at which dive started
#' @param d0.last If the depth bin that proceeded the first depth bin in 
#'   \code{depths}.  If the trajectory to be analyzed was started at the 
#'   surface, then set \code{c0.last=NULL}.
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#' 
#' @example examples/dsdive.ld.fixedstages.R
#' 
#' @importFrom stats dexp dbinom
#' 
#' @export
#' 
dsdive.ld.fixedstages = function(depths, durations, times, stages, beta, lambda, 
                                 depth.bins) {
  
  # extract dimensional information
  nt = length(times)
  
  # initialize log-density
  ld = 0
  
  # loop over transitions
  for(j in 1:max(nt-1, 1)) {
    
    # get transition parameters
    p = dsdive.tx.params(depth.bins = depth.bins, d0 = depths[j], 
                         s0 = stages[j], beta = beta, lambda = lambda)

    # account for duration
    # ld = ld + dexp(x = durations[j], rate = p$rate, log = TRUE)
    
    # only consider non-trivial depth bin transitions
    d0 = depths[j]
    df = depths[j+1]
    if(d0 != df) {
      if(length(p$labels) > 1) {
        went.down = df > d0
        ld = ld + dbinom(x = went.down, size = 1, prob = p$probs[p$labels > d0])
      }  
    }
  }
  
  
  #
  # account for durations
  #
  
  pace = durations / (2 * depth.bins[depths, 2])
  
  p1 = pace[stages==1]
  p2 = pace[stages==2]
  p3 = pace[stages==3]
  
  ld = ld + 
    dexp(x = sum(p1[is.finite(p1)&!is.na(p1)]), rate = lambda[1], log = TRUE) +
    dexp(x = sum(p2[is.finite(p2)&!is.na(p2)]), rate = lambda[2], log = TRUE) +
    dexp(x = sum(p3[is.finite(p3)&!is.na(p3)]), rate = lambda[3], log = TRUE)
    
  ld
}