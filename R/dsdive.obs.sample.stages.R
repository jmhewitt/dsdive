#' Sample stage transition times from full conditional posterior
#' 
#' Sampler uses a piecewise quadratic polynomial approximation to the log of the 
#' full conditional posterior.  Note that this sampler will sample both stage 
#' transition times at the same time.
#' 
#' @param depths Indices of observed depth bins
#' @param times Times at which each of \code{depths} was observed
#' @param t.stages Stage transition times for the dive; will be used to compute
#'   the dive stage for each observation
#' @param P.raw list of continuous time probability transition matrices, and 
#'   components.
#' @param T.range start and end times for the dive; used to control the support 
#'   of the stage transition time distributions
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param T1.prior.params \code{shape} and \code{rate} parameters for Gamma 
#'   prior on the descent-stage duration.
#' @param T2.prior.params \code{shape} and \code{rate} parameters for Gamma 
#'   prior on the bottom-stage duration.
#' @param max.width The stage transition times are updated with a piecewise 
#'   proposal distribution.  \code{max.width} controls the maximum width of the 
#'   intervals for the proposal distribution.  This is a tuning parameter that 
#'   controls the numerical stability of the proposal distribution, which is 
#'   sampled via inverse CDF techniques.
#' @param debug \code{TRUE} to return debugging objects, such as the proposal 
#'   density and log posterior
#'   
#' @example examples/dsdive.obs.sample.stages.R
#' 
#' @importFrom stats dgamma runif
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
  
  # can't allow stage 1->2 transition at surface
  inds.jumps = inds.jumps[depths[inds.jumps]!=1]
  
  # make sure inds.jumps exceed T.range
  inds.jumps = inds.jumps[times[inds.jumps] >= T.range[1]]
  
  # determine intervals and midpoints for log-quadratic sampling envelope
  breaks = refine.partition(
    breaks = sort(unique(c(T.range[1], times[inds.jumps], t.stages[2]))),
    max.width = max.width)
  anchors = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  
  # evaluate log-posterior at anchor and break points
  lp.anchors = lp(x = anchors, stage.ind = 1, t.stages = t.stages)
  lp.breaks = lp(x = breaks, stage.ind = 1, t.stages = t.stages)
  
  # build envelope
  envelope = envelope.approx(breaks = breaks, anchors = anchors, 
                             lp.breaks = lp.breaks, lp.anchors = lp.anchors)
  q1 = envelope.logquad(breaks = breaks, logf = envelope$logf, 
                        d.logf = envelope$d.logf, 
                        dd.logf.sup = envelope$dd.logf.sup, 
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
  
  # can't allow stage 2->3 transition at surface return
  surface.inds = inds.jumps[depths[inds.jumps]==1]
  if(length(surface.inds) > 0) {
    tmax = min(T.range[2], min(times[surface.inds]))
  } else {
    tmax = T.range[2]
  }
  inds.jumps = inds.jumps[depths[inds.jumps]!=1]
  
  # determine intervals and midpoints for log-quadratic sampling envelope
  breaks = refine.partition(
    breaks = sort(unique(c(t.stages[1], times[inds.jumps], tmax))),
    max.width = max.width)
  anchors = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  
  # evaluate log-posterior at anchor and break points
  lp.anchors = lp(x = anchors, stage.ind = 2, t.stages = t.stages)
  lp.breaks = lp(x = breaks, stage.ind = 2, t.stages = t.stages)
  
  
  # build envelope
  envelope = envelope.approx(breaks = breaks, anchors = anchors, 
                             lp.breaks = lp.breaks, lp.anchors = lp.anchors)
  q2 = envelope.logquad(
    breaks = breaks, logf = envelope$logf,
    d.logf = envelope$d.logf,
    dd.logf.sup = envelope$dd.logf.sup,
    anchors = anchors)

  # draw proposal
  prop = q2$rquad(n = 1)

  # compute metropolis ratio (as an independence sampler)
  if(t.stages[2] > T.range[2]) {
    lR = 0
  } else {
    lR = (lp(x = prop, stage.ind = 2, t.stages = t.stages) -
            q2$dquad(x = prop, log = TRUE)) -
      (lp(x = t.stages[2], stage.ind = 2, t.stages = t.stages) -
         q2$dquad(x = t.stages[2], log = TRUE))
  }
  
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