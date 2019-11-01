#' Likelihood for completely observed dive trajectories
#'
#' 
#' @example examples/ld.R
#' 
#' @export
#' 
#'
dsdive.ld = function(depths, durations, times, stages, beta, lambda, sub.tx,
                     surf.tx, depths.labels) {
  
  # extract dimensional information
  nt = length(times)
  num.depths = length(depths.labels) - 1
  
  # initialize log-density
  ld = 0
  
  # loop over transitions
  for(j in 1:(nt-1)) {
    
    # get transition parameters
    p = dsdive.tx.params(t0 = times[j], num.depths = num.depths, 
                         d0 = depths[j], s0 = stages[j], beta = beta, 
                         d0.last = ifelse(j==1, NULL, depths[j-1]), 
                         lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx)
    #
    # build likelihood
    #
    
    ld = ld +
      # duration
      dexp(x = durations[j], rate = p$rate, log = TRUE) + 
      # stage transition
      dbinom(x = stages[j] != stages[j+1], size = 1, prob = p$prob.stage, 
             log = TRUE)
    
    # add state transition
    if(length(p$labels) > 1) {
      ld = ld + log(p$probs[which(depths[j+1] == p$labels), stages[j+1]])
    }
    
  }
  
  ld
}