#' Likelihood for completely observed dive trajectories
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
#' @example examples/ld.R
#' 
#' @importFrom stats dexp dbinom
#' 
#' @export
#' 
dsdive.ld = function(depths, durations, times, stages, beta, lambda, sub.tx,
                     surf.tx, depth.bins, t0.dive, d0.last = NULL, t.stage2,
                     model) {
  
  # extract dimensional information
  nt = length(times)
  num.depths = nrow(depth.bins)
  
  # initialize log-density
  ld = 0
  
  # loop over transitions
  for(j in 1:max(nt-1,1)) {
    
    if(j>1) {
      d0.last = depths[j-1]
    }
    
    # get transition parameters
    p = dsdive.tx.params(t0 = times[j], depth.bins = depth.bins, 
                         d0 = depths[j], s0 = stages[j], beta = beta, 
                         d0.last = d0.last, 
                         lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx,
                         t0.dive = t0.dive, t.stage2 = t.stage2, model = model)
    #
    # build likelihood
    #
    
    # duration
    ld = ld + dexp(x = durations[j], rate = p$rate, log = TRUE)
    # stage transition
    if(length(stages) > 1) {
      ld = ld + dbinom(x = stages[j] != stages[j+1], size = 1, 
                       prob = p$prob.stage, log = TRUE)
      
      # add state transition
      if(length(p$labels) > 1) {
        ld = ld + log(p$probs[which(depths[j+1] == p$labels), stages[j+1]])
      }
    }
    
  }
  
  ld
}