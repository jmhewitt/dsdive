#' Simulate dive trajectories across discrete depth bins
#'
#' The method will simulate dive trajectories from initial conditions until the 
#' trajectory is observable at \code{tf}, or a maximum number of transitions 
#' has been exceeded.  The dive simulation is bridged, so the trajectory will
#' also stop diving after returning to the surface.
#' 
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param d0 the depth bin at which transition parameters should be computed
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
#' @param t0 time at which transition parameters should be computed
#' @param tf time at which sampling should end after
#' @param s0 dive stage in which forward simulation begins
#' 
#' @return A \code{dsdive} object, which is a \code{list} with the following 
#'   vectors:
#'   \describe{
#'     \item{depths}{Record of which depth bin indices the trajectory visited}
#'     \item{durations}{Record of amount of time spent in each depth bin}
#'     \item{times}{The time at which each depth bin was entered}
#'     \item{stages}{The stage at which each depth bin was entered}
#'   }
#' 
#' @example examples/dsdive.impute.R
#' 
#' @export
#'
dsdive.impute = function(depth.bins, depths, times, s0, beta, lambda, sub.tx, 
                         surf.tx, inflation.factor.lambda = 1.1, 
                         verbose = FALSE, precompute.bridges = TRUE) {
  
  # extract dimensions
  nt = length(times)
  
  # initialize imputed trajectory
  res = list(
    depths = depths[1],
    stages = s0,
    times = times[1],
    durations = c(),
    ld = 0
  )
  
  # impute segments
  for(i in 1:(nt-1)) {
    
    # extract last depth bin index
    if(i==1) { 
      d0.last = NULL 
    } else { 
      d0.last = res$depths[length(res$depths) - 1]
    } 
    
    # impute trajectory segment
    br = dsdive.bridgesample(depth.bins = depth.bins, d0 = depths[i], 
                             d0.last = d0.last, df = depths[i+1], beta = beta, 
                             lambda = lambda, sub.tx = sub.tx, 
                             surf.tx = surf.tx, t0 = times[i], tf = times[i+1], 
                             s0 = res$stages[length(res$stages)], 
                             inflation.factor.lambda = inflation.factor.lambda,
                             verbose = verbose, 
                             precompute.bridges = precompute.bridges)
    
    # ensure initial imputed duration yields continuously observable trajectory
    if(i > 1) {
      br$durations[1] = 
        br$durations[1] + times[i] - res$times[length(res$times)]
    }
    
    # merge segments
    res$depths = c(res$depths, br$depths[-1])
    res$stages = c(res$stages, br$stages[-1])
    res$times = c(res$times, br$t[-1])
    res$durations = c(res$durations, br$durations)
    
    # update log-density for proposed trajectory
    res$ld = res$ld + br$ld
    
  }
  
  class(res) = 'dsdive'
  
  res
}