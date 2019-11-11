#' Use bridged sampling to impute a complete dive trajectory consistent with observations
#'
#' 
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param depths record of depth bins the trajectory should visit
#' @param times times at which the depth bins should be visited
#' @param s0 dive stage at which the trajectory should be started from
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
#' @param inflation.factor.lambda In order to facilitate bridged transitions, 
#'   the transition rate of the overall process must be inflated to allow the 
#'   possibility of self-transitions.  Self-transitions allow bridged paths to 
#'   dynamically modify the total number of transitions between observed values
#'   so that a valid path between observations is always possible.  The 
#'   \code{inflation.factor.lambda} parameter implicitly controls the number of 
#'   self-transitions that will occur.  Larger values will create more 
#'   self-transitions.
#' @param verbose If \code{TRUE}, then the sampler's progress will be printed 
#'   during sampling.
#' @param precompute.bridges If \code{TRUE}, then the bridged transition 
#'   matrices will be precomputed.  Enabling this option will increase the 
#'   memory overhead of the method, but will reduce its runtime.
#' @param t0.dive Time at which dive started
#' 
#' @example examples/dsdive.impute.R
#' 
#' @export
#'
dsdive.impute = function(depth.bins, depths, times, s0, beta, lambda, sub.tx, 
                         surf.tx, inflation.factor.lambda = 1.1, 
                         verbose = FALSE, precompute.bridges = TRUE, t0.dive) {
  
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
  
  # pre-compute the common max transition rate
  lambda.max = max(outer(lambda, 2 * depth.bins[,2], '/'))
  
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
                             precompute.bridges = precompute.bridges,
                             lambda.max = lambda.max, t0.dive = t0.dive)
    
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