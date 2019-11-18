#' Use bridged sampling to impute a complete dive trajectory consistent with observations
#'
#' The sampling method is designed to sample many trajectories simultaneously, 
#' so has an extra level of approximation in the proposal distributions.
#' 
#' @param M the number of trajectories to sample
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
#' @param resample Resample particles at each step if \code{TRUE}.
#' @param trajectory.conditional If not \code{NULL}, then 
#'   \code{trajectory.conditional} must be the dive information for a 
#'   completely observed \code{dsdive} object.  The first entries in the 
#'   initialization vectors \code{d0}, \code{d0.last}, \code{s0} must be 
#'   associated with the trajectory observed in \code{trajectory.conditional}.
#'   Providing a non \code{NULL} value for \code{trajectory.conditional} will 
#'   cause \code{dsdive.fastbridge} to simulate one fewer trajectories, as 
#'   the value of \code{trajectory.conditional} will be returned.  The 
#'   importance of this function argument is that \code{dsdive.fastbridge}
#'   will evaluate the proposal density for \code{trajectory.conditional}.
#' 
#' @example examples/dsdive.impute.R
#' 
#' @export
#'
dsdive.fastimpute = function(M, depth.bins, depths, times, s0, beta, lambda, 
                             sub.tx, surf.tx, inflation.factor.lambda = 1.1, 
                             verbose = FALSE, precompute.bridges = TRUE, 
                             t0.dive, resample = TRUE, 
                             trajectory.conditional = NULL) {
  
  # set flag for conditional simulation
  cond.sim = !is.null(trajectory.conditional)
  if(cond.sim) {
    resample = TRUE
  }
  
  # extract dimensions
  nt = length(times)
  
  # initialize imputed trajectories
  res.single = list(
    depths = depths[1],
    stages = s0,
    times = times[1],
    durations = c(),
    ld = 0,
    ld.true = 0,
    w = 0
  )
  
  class(res.single) = 'dsdive'
  
  # initialize results
  res = replicate(M, list(res.single))
  
  # pre-compute the common max transition rate
  lambda.max = max(outer(lambda, 2 * depth.bins[,2], '/'))
  
  # impute segments
  for(i in 1:(nt-1)) {
    
    # extract last depth bin indices
    if(i==1) { 
      d0.last = NULL
    } else { 
      d0.last = sapply(res, function(r) r$depths[length(r$depths) - 1])
    } 
    
    # extract current state indices
    s0.current = sapply(res, function(r) r$stages[length(r$stages)])
    
    # set conditional simulation parameters
    if(cond.sim) {
      segment.conditional = dsdive.extract(
        depths = trajectory.conditional$depths, 
        times = trajectory.conditional$times, 
        stages = trajectory.conditional$stages, 
        durations = trajectory.conditional$durations, 
        t0 = times[i], tf = times[i+1])
    } else {
      segment.conditional = NULL
    }
    
    # impute trajectory segment
    br = dsdive.fastbridge(M = M, depth.bins = depth.bins, d0 = depths[i], 
                           d0.last = d0.last, df = depths[i+1], beta = beta, 
                           lambda = lambda, sub.tx = sub.tx, 
                           surf.tx = surf.tx, t0 = times[i], tf = times[i+1], 
                           s0 = s0.current,
                           inflation.factor.lambda = inflation.factor.lambda,
                           verbose = verbose, 
                           precompute.bridges = precompute.bridges,
                           lambda.max = lambda.max, t0.dive = t0.dive, 
                           trajectory.conditional = segment.conditional)
    
    # ensure initial imputed duration yields continuously observable trajectory
    if(i > 1) {
      for(j in 1:M) {
        if(!(cond.sim & j==1)) {
          br[[j]]$durations[1] = br[[j]]$durations[1] + times[i] - 
            res[[j]]$times[length(res[[j]]$times)]
        }
      }
    }
  
    for(j in 1:M) {
      
      # merge segments
      res[[j]]$depths = c(res[[j]]$depths, br[[j]]$depths[-1])
      res[[j]]$stages = c(res[[j]]$stages, br[[j]]$stages[-1])
      res[[j]]$times = c(res[[j]]$times, br[[j]]$t[-1])
      if(length(br[[j]]$depths) > 1) {
        res[[j]]$durations = c(res[[j]]$durations, br[[j]]$durations)
      }
      
      # update log-density for proposal
      res[[j]]$ld = res[[j]]$ld + br[[j]]$ld
      
      # compute sampling log-weight
      res[[j]]$w = br[[j]]$ld.true - br[[j]]$ld
    }
    
    # resample particles
    if(resample) {
      # get raw weights
      W = exp(sapply(res, function(r) r$w))
      # correct for edge cases
      W[is.infinite(W)] = sign(W[is.infinite(W)])
      W[W==-1] = 0
      # standardize weights
      W = W / sum(W)
      if(cond.sim) {
        res = c(res[1], res[sample(x = M, size = M - 1, replace = TRUE, 
                                   prob = W)])
      } else{
        res = res[sample(x = M, size = M, replace = TRUE, prob = W)]
      }
    }
    
  }
  
  # compute log-densities under true model
  for(j in 1:M) {
    rj = res[[j]]
    res[[j]]$ld.true = dsdive.ld(depths = rj$depths, durations = rj$durations, 
                                 times = rj$times, stages = rj$stages,
                                 beta = beta, lambda = lambda, sub.tx = sub.tx, 
                                 surf.tx = surf.tx, depth.bins = depth.bins, 
                                 t0.dive = t0.dive)
  }
  
  res
}