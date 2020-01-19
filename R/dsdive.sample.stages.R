#' Sample stage transition times from full conditional posterior
#' 
#' @param prior.t list of two functions that will allow the prior density for 
#'   stage transition times to be called
#' 
#' @export
#' 
dsdive.sample.stages = function(depths, durations, times, t.stages, 
                                beta, lambda, depth.bins, T1.prior, T2.prior) {
  
  # get sequence of downward transitions
  tx.down = diff(depths) == 1
  
  # get indices of transitions
  tx.inds = 1:length(tx.down)
  
  # get widths of all depth bins
  widths = 2 * depth.bins[depths[tx.inds], 2]
  
  # compute dive stages, given conditional information
  stages = findInterval(times[tx.inds], t.stages) + 1
  
  #
  # build conditional likelihood for stage 1->2 transition time
  #
  
  s12 = stages<3
  
  # density of each observation under stage 1/2 parameters
  ld.1 = dbinom(x = tx.down[s12], size = 1, prob = beta[1], log = TRUE) + 
    dexp(x = durations[s12], rate = lambda[1] / widths[s12], log = TRUE)
  ld.2 = -0.6931472 + 
    dexp(x = durations[s12], rate = lambda[2] / widths[s12], log = TRUE)
  
  s3.start = min(which(!s12))
  if(is.infinite(s3.start)) {
    s3.start = length(ld.2) - 1
  }
  
  # determine support of stage 1->2 transition
  support.tx12 = 2:max(s3.start - 1, 2)
  
  # determine log-posterior across support
  lp = T1.prior
  d = numeric(length(support.tx12))
  nd = length(d)
  nl2 = length(ld.2)
  for(i in 1:nd) {
    b = support.tx12[i]
    d[i] = sum(ld.1[1:(b-1)]) + sum(ld.2[b:nl2]) + lp((times[b] - times[1])/60)
  }
  
  # sample new stage 1->2 transition index/time
  tx.12.ind = sample.gumbeltrick(log.p = d)
  
  # build new stage vector
  stages = stagevec(length.out = length(stages), 
                    breaks = c(tx.12.ind, s3.start))
  
  #
  # build conditional likelihood for stage 2->3 transition time
  #
  
  s23 = stages>1
  
  # density of each observation under stage 2/3 parameters
  ld.2 = -0.6931472 + 
    dexp(x = durations[s23], rate = lambda[2] / widths[s23], log = TRUE)
  ld.3 = dbinom(x = tx.down[s23], size = 1, prob = beta[2], log = TRUE) + 
    dexp(x = durations[s23], rate = lambda[3] / widths[s23], log = TRUE)
  
  # determine support of stage 1->2 transition
  nt = length(times)
  if(tx.12.ind+1 < nt) {
    support.tx23 = (tx.12.ind+1):nt
  } else {
    support.tx23 = nt
  }
  
  # determine log-posterior across support
  t.stage2 = times[tx.12.ind]
  lp = T2.prior
  d = numeric(length(support.tx23))
  nd = length(d)
  nl3 = length(ld.3)
  for(i in 1:nd) {
    d[i] = sum(ld.2[1:(i-1)]) + sum(ld.3[(i+1):nl3]) + 
      lp((times[support.tx23[i]] - t.stage2)/60)
  }
  
  # sample new stage 2->3 transition index/time
  tx.23.ind = support.tx23[sample.gumbeltrick(log.p = d)]
  
  # build new stage vector
  stages = c(stagevec(length.out = length(times), 
                    breaks = c(tx.12.ind, tx.23.ind)))
  
  # package results
  list(
    stages = stages,
    t.stages = times[c(tx.12.ind, tx.23.ind)]
  )
}