#' Sample stage transition times from full conditional posterior
#' 
#' @param prior.t list of two functions that will allow the prior density for 
#'   stage transition times to be called
#' 
#' @example examples/dsdive.sample.stages.R
#' 
#' @export
#' 
dsdive.sample.stages = function(depths, durations, times, t.stages, 
                                beta, lambda, depth.bins, T1.prior, T2.prior) {
  
  n = length(depths)
  
  #
  # sample stage 1->2 transition, conditional on stage 3 time
  #
  
  # compute dive stages, given conditional information
  stages = findInterval(times, t.stages) + 1
  
  # get index of stage 3 time
  s3.start = min(which(stages == 3))
  if(is.infinite(s3.start)) {
    s3.start = n-1
  }
  
  # determine support of stage 1->2 transition
  s2.start.support = 2:max(s3.start-1, 2)
  ns2 = length(s2.start.support)
  
  # density of each transition index under stage 1 and stage 2 parameters
  ld1 = numeric(length = ns2)
  ld2 = numeric(length = ns2)
  for(i in 1:ns2) {
    
    # map back to a depth bin in the dive
    depth.ind = s2.start.support[i]
    d0 = depths[depth.ind]
    tau = durations[depth.ind]
    
    # stage 1 and stage 2 parameters for depth bin
    
    p1 = dsdive.tx.params(depth.bins = depth.bins, d0 = d0, s0 = 1, 
                          beta = beta, lambda = lambda)
    
    p2 = dsdive.tx.params(depth.bins = depth.bins, d0 = d0, s0 = 2, 
                          beta = beta, lambda = lambda)
    
    # duration component of likelihood
    if(is.finite(tau)) {
      ld1[i] = dexp(x = tau, rate = p1$rate, log = TRUE)
      ld2[i] = dexp(x = tau, rate = p2$rate, log = TRUE)
    }
    
    # depth bin transition component of likelihood
    if(depth.ind + 1 <= n) {
      df = depths[depth.ind + 1]
      success1 = p1$labels[1] == df
      success2 = p2$labels[1] == df
      ld1[i] = ld1[i] + dbinom(x = success1, size = 1, prob = p1$probs[1], 
                               log = TRUE)
      ld2[i] = ld2[i] + dbinom(x = success2, size = 1, prob = p2$probs[1], 
                               log = TRUE)
    }
    
  }
  
  # compute log-posterior across support
  d = numeric(length = ns2)
  for(i in 1:ns2) {
    d[i] = ifelse(i>1, sum(ld1[1:(i-1)]), 0) + sum(ld2[i:ns2]) +
      T1.prior((times[s2.start.support[i]] - times[1])/60)
  }
  
  # sample new time at which stage 2 begins
  s2.start = s2.start.support[sample.gumbeltrick(log.p = d)]
  
  # build new stage vector
  stages = stagevec(length.out = length(depths), breaks = c(s2.start, s3.start))
  
  
  #
  # sample stage 2->3 transition, conditional on stage 2 time
  #

  # determine support of stage 2->3 transition
  s3.start.support = (s2.start+1):(n-1)
  ns3 = length(s3.start.support)
  
  # density of each transition index under stage 2 and stage 3 parameters
  ld2 = numeric(length = ns3)
  ld3 = numeric(length = ns3)
  for(i in 1:ns3) {
    
    # map back to a depth bin in the dive
    depth.ind = s3.start.support[i]
    d0 = depths[depth.ind]
    tau = durations[depth.ind]
    
    # stage 2 and stage 3 parameters for depth bin
    
    p2 = dsdive.tx.params(depth.bins = depth.bins, d0 = d0, s0 = 2, 
                          beta = beta, lambda = lambda)
    
    p3 = dsdive.tx.params(depth.bins = depth.bins, d0 = d0, s0 = 3, 
                          beta = beta, lambda = lambda)
    
    # duration component of likelihood
    if(is.finite(tau)) {
      ld2[i] = dexp(x = tau, rate = p2$rate, log = TRUE)
      ld3[i] = dexp(x = tau, rate = p3$rate, log = TRUE)
    }
    
    # depth bin transition component of likelihood
    if(depth.ind + 1 <= n) {
      df = depths[depth.ind + 1]
      success2 = p2$labels[1] == df
      success3 = p3$labels[1] == df
      ld2[i] = ld2[i] + dbinom(x = success2, size = 1, prob = p2$probs[1], 
                               log = TRUE)
      ld3[i] = ld3[i] + dbinom(x = success3, size = 1, prob = p3$probs[1], 
                               log = TRUE)
    }
    
  }
  
  # compute log-posterior across support
  d = numeric(length = ns3)
  for(i in 1:ns3) {
    d[i] = ifelse(i>1, sum(ld2[1:(i-1)]), 0) + sum(ld3[i:ns3]) + 
      T2.prior((times[s3.start.support[i]] - times[s2.start])/60)
  }
  
  # sample new time at which stage 3 begins
  s3.start = s3.start.support[sample.gumbeltrick(log.p = d)]
  
  # build new stage vector
  stages = stagevec(length.out = length(depths), breaks = c(s2.start, s3.start))
  
  
  # package results
  list(
    stages = stages,
    t.stages = times[c(s2.start, s3.start)]
  )
}