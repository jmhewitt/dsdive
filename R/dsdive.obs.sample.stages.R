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
dsdive.obs.sample.stages = function(depths, times, t.stages, P.raw, 
                                    T.range, depth.bins, T1.prior.params, 
                                    T2.prior.params, max.width, debug = FALSE) {
  
  # package into dsobs object, for likelihood evaluation
  dsobs = list(depths = depths, times = times)
  class(dsobs) = 'dsobs'
  
  
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
    
    # log density associated with stage transition time
    ld = dsdive.obsld(dsobs.list = dsobs, t.stages.list = tstages, 
                      P.raw = P.raw, s0 = stage.ind, sf = stage.ind + 1)
    
    # extract the parameters for the prior distribution, and convert the stage 
    # transition time into a duration
    if(stage.ind==1) { 
      prior.params = T1.prior.params 
      x.dur = x - T.range[1]
    } else if(stage.ind==2) { 
      prior.params = T2.prior.params 
      x.dur = x - t.stages[1]
    }
    
    # return log-posterior by adding log prior
    ld + dgamma(x = x.dur, shape = prior.params[1], rate = prior.params[2], 
                log = TRUE)
  })}
  
  
  #
  # sample stage 1->2 transition time
  #
  
  # identify timepoints where there are jump-discontinuities in log posterior
  inds.jumps = which(times <= t.stages[2])
  
  # determine intervals and midpoints for log-quadratic sampling envelope
  breaks = refine.partition(
    breaks = sort(unique(c(T.range[1], times[inds.jumps], t.stages[2]))),
    max.width = max.width)
  anchors = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  
  # evaluate log-posterior at anchor points
  anchors.extra = c(anchors, breaks[length(breaks)])
  lp.eval = lp(x = anchors.extra, stage.ind = 1, t.stages = t.stages)
  
  # build envelope
  q1 = envelope.logquad(breaks = breaks, logf = lp.eval[1:length(anchors)], 
                        d.logf = diff(lp.eval)/diff(anchors.extra), 
                        dd.logf.sup = rep(0, length(anchors)), 
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
    breaks = sort(unique(c(t.stages[1], times[inds.jumps], T.range[2]))),
    max.width = max.width)
  anchors = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  
  # evaluate log-posterior at anchor points
  anchors.extra = c(anchors, breaks[length(breaks)])
  lp.eval = lp(x = anchors.extra, stage.ind = 2, t.stages = t.stages)

  # build envelope
  q2 = envelope.logquad(
    breaks = breaks, logf = lp.eval[1:length(anchors)],
    d.logf = diff(lp.eval)/diff(anchors.extra),
    dd.logf.sup = rep(0, length(anchors)),
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
    res$debug = list(lp = lp, q1 = q1, q2 = q2)
  }
  
  res
}