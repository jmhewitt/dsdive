#' Sample stage transition times from full conditional posterior
#' 
#' We assume the prior distributions for the stage transition times are Gamma 
#' distributions with shape/rate parameterization.
#' 
#' @param max.width for numerical stability of the sampling envelopes
#' 
#' @example examples/dsdive.sample.stages.R
#' 
#' @export
#' 
dsdive.sample.stages = function(depths, times, t.stages, beta, lambda, 
                                depth.bins, T1.prior.params, T2.prior.params,
                                max.width, debug = FALSE) {
  
  #
  # components for constructing proposal densities for full conditionals
  # 
  
  lp = function(x, stage.ind, t.stages) { sapply(x, function(x) {
    # log-posterior for data given a stage transition time
    # 
    # Parameters: 
    #  x - stage transition time
    #  stage.ind - 1 if x represents the 1->2 stage transition time; 2 otherwise
    #  t.stages - current values of stage transition times
    
    tstages = t.stages
    tstages[stage.ind] = x
    
    # split depth bin durations according to the stage transition times
    augmented = dsdive.augment.trajectory(depths = depths, times = times, 
                                          t.stages = tstages)
    
    # log density associated with stage transition time
    ld = dsdive.ld.fixedstages(
      depths = augmented$depths, durations = augmented$durations, 
      times = augmented$times, stages = augmented$stages, beta = beta, 
      lambda = lambda, depth.bins = depth.bins)
    
    # extract the parameters for the prior distribution, and convert the stage 
    # transition time into a duration
    if(stage.ind==1) { 
      prior.params = T1.prior.params 
      x.dur = x - times[1]
    } else if(stage.ind==2) { 
      prior.params = T2.prior.params 
      x.dur = x - t.stages[1]
    }
    
    # return log-posterior by adding log prior
    ld + dgamma(x = x.dur, shape = prior.params[1], rate = prior.params[2], 
                log = TRUE)
  })}
  
  dlp = function(x, stage.ind, t.stages) { sapply(x, function(x) {
    # derivative of log-posterior given a stage transition time
    # 
    # Parameters:
    #  x - stage transition time
    #  stage.ind - 1 if x represents the 1->2 stage transition time; 2 otherwise
    #  t.stages - current values of stage transition times
    
    # get depth bin at time x
    bin = depths[findInterval(x, times)]
    
    # extract the parameters for the prior distribution, and convert the stage 
    # transition time into a duration
    if(stage.ind==1) { 
      prior.params = T1.prior.params 
      x.dur = x - times[1]
    } else if(stage.ind==2) { 
      prior.params = T2.prior.params 
      x.dur = x - t.stages[1]
    }
    
    diff(lambda[stage.ind + 0:1]) / (2*depth.bins[bin, 2]) - prior.params[2] + 
      (prior.params[1] - 1) / x.dur
  })}
  
  ddlp.sup = function(breaks, stage.ind, t.stages) {
    # suprema of second derivative for log-posterior across intervals defined 
    # via breaks
    #
    # Parameters:
    #  breaks - endpoints of intervals
    #  stage.ind - 1 if x represents the 1->2 stage transition time; 2 otherwise
    #  t.stages - current values of stage transition times
    
    # extract the parameters for the prior distribution, and convert the 
    # interval times into durations
    if(stage.ind==1) { 
      prior.params = T1.prior.params 
      x.dur = breaks - times[1]
    } else if(stage.ind==2) { 
      prior.params = T2.prior.params 
      x.dur = breaks - t.stages[1]
    }
    
    n = length(breaks)
    if(prior.params[1] == 0) {
      r = rep(0, n-1)
    } else if(prior.params[1] > 0) {
      r = (1 - prior.params[1]) / (x.dur[-1])^2
    } else if(prior.params[1] < 0) {
      r = (1 - prior.params[1]) / (x.dur[1:(n-1)])^2
    }
    
    # zap small curvatures to increase numerical stability of envelope sampler.
    # Notes: 
    #  1) If the curvature was positive, then we might not have a true envelope
    #  2) For reasonable prior parameterizations, we will have neg. curvature
    r[abs(r) < 1e-5] = 0
    
    r
  }
  
  
  #
  # sample stage 1->2 transition time
  #
  
  # identify timepoints where there are jump-discontinuities in log posterior
  inds.jumps = which(times <= t.stages[2])
  
  # determine intervals and midpoints for log-quadratic sampling envelope
  breaks = refine.partition(
    breaks = sort(unique(c(times[1], times[inds.jumps], t.stages[2]))),
    max.width = max.width)
  anchors = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  
  # build envelope
  q1 = envelope.logquad(
    breaks = breaks, logf = lp(x = anchors, stage.ind = 1, t.stages = t.stages), 
    d.logf = dlp(x = anchors, stage.ind = 1, t.stages = t.stages), 
    dd.logf.sup = ddlp.sup(breaks = breaks, stage.ind = 1, t.stages = t.stages), 
    anchors = anchors)
  
  # draw proposal
  prop = q1$rquad(n = 1)
  
  # compute metropolis ratio (as an independence sampler)
  lR = (lp(x = prop, stage.ind = 1, t.stages = t.stages) - 
          q1$dquad(x = prop, log = TRUE)) - 
       (lp(x = t.stages[1], stage.ind = 1, t.stages = t.stages) - 
          q1$dquad(x = t.stages[1], log = TRUE))
  
  # accept/reject
  if(log(runif(1)) <= lR) {
    t.stages[1] = prop
  }
    
  
  #
  # sample stage 2->3 transition time
  #
  
  # identify timepoints where there are jump-discontinuities in log posterior
  inds.jumps = which(times >= t.stages[1])
  
  # determine intervals and midpoints for log-quadratic sampling envelope
  breaks = refine.partition(
    breaks = sort(unique(c(t.stages[1], times[inds.jumps], 
                           times[length(times)]))),
    max.width = max.width)
  anchors = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  
  # build envelope
  q2 = envelope.logquad(
    breaks = breaks, logf = lp(x = anchors, stage.ind = 2, t.stages = t.stages), 
    d.logf = dlp(x = anchors, stage.ind = 2, t.stages = t.stages), 
    dd.logf.sup = ddlp.sup(breaks = breaks, stage.ind = 2, t.stages = t.stages), 
    anchors = anchors)
  
  # draw proposal
  prop = q2$rquad(n = 1)
  
  # compute metropolis ratio (as an independence sampler)
  lR = (lp(x = prop, stage.ind = 2, t.stages = t.stages) - 
          q2$dquad(x = prop, log = TRUE)) - 
       (lp(x = t.stages[2], stage.ind = 2, t.stages = t.stages) - 
          q2$dquad(x = t.stages[2], log = TRUE))
  
  # accept/reject
  if(log(runif(1)) <= lR) {
    t.stages[2] = prop
  }
  
  
  #
  # Package results
  #
  
  res = list(t.stages = t.stages)
  
  if(debug == TRUE) {
    res$debug = list(lp = lp, dlp = dlp, ddlp.sup = ddlp.sup, q1 = q1, q2 = q2)
  }
  
  res
}