#' Impute a complete dive trajectory from partial observations
#' 
#' Uses properties of homogeneous Continuous time Markov Chains (CTMCs) to 
#' impute trajectory segments between observations.  The sampler begins by 
#' sampling the number of transitions between observations \eqn{N}, then 
#' samples a length \code{N} path that connects the trajectory segment at its 
#' start and end states.
#' 
#' @param depths Depth bin indices visited
#' @param times Times at which each of \code{depths} was visited
#' @param t.stages Stage transition times for the dive; will be used to compute
#'   the dive stage for each observation
#' @param rate.unif uniformization rate, for standardizing transition
#'   rates between states
#' @param P.raw list of continuous time probability transition matrices, and 
#'  components.
#' @param P.tx list of discrete time probability transition matrices
#' @param n.bins number of rows in the \code{depths} matrix
#' @param max.tx maximum number of transitions between observations that will 
#'   be allowed during imputation
#' 
#' @example examples/dsdive.impute.R
#' 
#' @importFrom stats runif
#' 
#' @export
#'
dsdive.impute = function(depths, times, t.stages, rate.unif, P.raw, 
                         P.tx, n.bins = nrow(depths), max.tx) {
  
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
    
    # number of transitions to impute
    n = dsdive.impute.sample_n(
      d0 = d0, df = df, s0 = s0, sf = sf, t0 = t0, tf = tf, 
      t.stages = t.stages, rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, 
      n.bins = n.bins, max.tx = max.tx)
    
    
    # assemble single-step transition matrices
    s.range = s0:sf
    B = do.call(c, lapply(1:length(s.range), function(i) {
      s = s.range[i]
      if(n[i] > 0) {
        lapply(1:n[i], function(j) P.tx[[s]])
      } else {
        list()
      }
    }))
    
    # assemble likelihood for transitions
    L = matrix(0, nrow = n.bins, ncol = sum(n) + 1)
    L[d0,1] = 1                    # initial location is fixed
    L[df,ncol(L)] = 1              # final location is fixed
    L[,-c(1,ncol(L))] = 1/n.bins   # all other transitions are free
    
    # impute via backward sampling
    depths.raw[[i]] = ffbs(B = B, L = L)
    
    # uniformly sample transition times
    times.raw[[i]] = t0 + c(0, sort((tf-t0) * runif(n = sum(n))))
  
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