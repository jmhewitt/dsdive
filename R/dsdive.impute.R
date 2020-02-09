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
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#' @param ld.compute \code{TRUE} to compute likelihood values as well.  This 
#'   is required if resampling or conditional trajectory imputation is used.
#' @param t.stages the times at which stage transitions occur
#' 
#' @example examples/dsdive.impute.R
#' @export
#'
dsdive.impute = function(depths, times, t.stages, rate.unif, P.raw, P.tx, 
                         n.bins, max.tx) {
  
  # number of observations
  nt = length(depths)
  
  # associate stages with observations
  stages = findInterval(times, t.stages) + 1
  
  # initialize containers for segments
  depths.raw = vector('list', nt-1)
  times.raw = vector('list', nt-1)
  
  # impute segments
  for(i in 1:(nt-1)) {
    
    # extract start and end information for segment
    d0 = depths[i]
    df = depths[i+1]
    s0 = stages[i]
    sf = stages[i+1]
    t0 = times[i]
    tf = times[i+1]
    
    if(s0==sf) {
      
      # number of transitions to impute
      n = dsdive.impute.sample_n(
        n0 = NULL, d0 = d0, df = df, s0 = s0, sf = sf, t0 = t0, tf = tf, 
        t.stages = t.stages, rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, 
        ff.s0 = NULL, bf.sf = NULL, n.bins = n.bins, max.tx = max.tx, 
        filters.out = TRUE)
      
      # assemble single-step transition matrices
      B = lapply(1:n$n, function(j) P.tx[[s0]])
      
      # assemble likelihood for transitions
      L = matrix(0, nrow = n.bins, ncol = n$n + 1)
      L[d0,1] = 1                    # initial location is fixed
      L[df,ncol(L)] = 1              # final location is fixed
      L[,-c(1,ncol(L))] = 1/n.bins   # all other transitions are free

      # impute via backward sampling
      depths.raw[[i]] = bs(a = n$ff.s0, B = B, L = L)
      
      # uniformly sample transition times
      times.raw[[i]] = t0 + c(0, sort((tf-t0) * runif(n = n$n)))
      
    } else {
      
      # number of stage s0 transitions to impute
      n0 = dsdive.impute.sample_n(
        n0 = NULL, d0 = d0, df = df, s0 = s0, sf = sf, t0 = t0, tf = tf, 
        t.stages = t.stages, rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, 
        ff.s0 = NULL, bf.sf = NULL, n.bins = n.bins, max.tx = max.tx, 
        filters.out = TRUE)
      
      # number of stage sf transitions to impute
      n1 = dsdive.impute.sample_n(
        n0 = n0$n, d0 = d0, df = df, s0 = s0, sf = sf, t0 = t0, tf = tf, 
        t.stages = t.stages, rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, 
        ff.s0 = n0$ff.s0, bf.sf = NULL, n.bins = n.bins, max.tx = max.tx, 
        filters.out = FALSE)
      
      # assemble single-step transition matrices
      B = list()
      if(n0$n > 0) { B = c(B, lapply(1:n0$n, function(j) P.tx[[s0]])) }
      if(n1 > 0) { B = c(B, lapply(1:n1, function(j) P.tx[[sf]])) }
      
      # assemble likelihood for transitions
      L = matrix(0, nrow = n.bins, ncol = n0$n + n1 + 1)
      L[d0,1] = 1                    # initial location is fixed
      L[df,ncol(L)] = 1              # final location is fixed
      L[,-c(1,ncol(L))] = 1/n.bins   # all other transitions are free
      
      # impute via backward sampling
      depths.raw[[i]] = ffbs(B = B, L = L)
      
      # uniformly sample transition times
      times.raw[[i]] = t0 + c(0, sort((tf-t0) * runif(n = n0$n + n1)))
      
    }
    
  }
  
  # assemble uniformized dive
  depths.full = do.call(c, depths.raw)
  times.full = do.call(c, times.raw)
  stages.full = findInterval(times.full, t.stages) + 1
    
  #
  # package complete dive
  #
  
  # identify depth bin or stage transition times
  tx.real = c(TRUE, diff(depths.full) != 0) | c(TRUE, diff(stages.full) == 1)
  
  res = list(
    depths = depths.full[tx.real],
    stages = stages.full[tx.real],
    times = times.full[tx.real]
  )
  
  res$durations = c(diff(res$times), 
                    times[length(times)] - res$times[length(res$times)])
  
  class(res) = 'dsdive'
  
  res
}