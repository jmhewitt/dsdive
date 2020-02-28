#' Sample stage transition times from full conditional posterior
#' 
#' We assume the prior distributions for the stage transition times are Gamma 
#' distributions with shape/rate parameterization.
#' 
#' @example examples/dsdive.obs.sample.offsets.R
#' 
#' 
dsdive.obs.sample.offsets = function(dsobs.aligned, dsobs.unaligned, offset, 
                                     t.stages, P.raw, depth.bins, tstep,
                                     max.width, t0.prior.params, 
                                     debug = FALSE) {
  
  tstep2 = 2*tstep
  
  #
  # components for constructing proposal densities for full conditionals
  # 
  
  align.obs = function(eps) {
    # construct a dsobs object with realigned observation times
    #
    # Parameters:
    #   eps - offset
    
    # realign observation times with beginning of dive
    times = dsobs.unaligned$times - eps
    
    # remove observations from before the dive begins
    nonneg.times = times >= 0
    times = times[nonneg.times]
    depths = dsobs.unaligned$depths[nonneg.times]
    
    # dives start at the surface; add this if necessary
    if(times[1] > 0) {
      times = c(0, times)
      depths = c(1, depths)
    }
    
    # package dive
    aligned = list(
      times = times,
      depths = depths,
      stages = findInterval(times, t.stages) + 1
    )
    class(aligned) = 'dsobs'
    
    aligned
  }
  
  
  lp = function(eps) { sapply(eps, function(eps) {
    # log-posterior for data given a dive offset
    # 
    # Parameters: 
    #  eps - offset
    
    # construct dsobs object with realigned observation times
    aligned = align.obs(eps)
    
    # log density associated with dive offset
    ld = dsdive.obsld(dsobs.list = aligned, t.stages.list = t.stages, 
                      P.raw = P.raw, s0 = 1, sf = 3)
    
    # return log-posterior (ignoring jacobian)
    ld + dbeta(x = (eps + tstep)/tstep2, 
               shape1 = t0.prior.params[1], shape2 = t0.prior.params[2], 
               log = TRUE)
  })}
  
  
  #
  # sample dive offset
  #
  
  # determine intervals and midpoints for log-quadratic sampling envelope
  breaks = refine.partition(
    breaks = c(-tstep, tstep),
    max.width = max.width)
  anchors = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  
  # evaluate log-posterior at all points; account for Beta prior boundary issues
  tol.boundary = .Machine$double.eps * tstep2
  end.inds = c(1,length(breaks))
  lp.breaks = lp(eps = c(breaks[1] + tol.boundary, 
                         breaks[-end.inds],
                         breaks[length(breaks)] - tol.boundary))
  lp.anchors = lp(eps = anchors)
  
  # approximate slopes, avoiding issues at endpoints
  d.logf = c(
    (lp.breaks[2] - lp.anchors[1])/(breaks[2]-anchors[1]),
    diff(lp.breaks[-end.inds])/diff(breaks[-end.inds]),
    (lp.anchors[length(lp.anchors)] - lp.breaks[length(lp.breaks)-1])/(
      anchors[length(anchors)] - breaks[length(breaks)-1]
    )
  )
  
  # approximate curvatures by assuming local quadratic fits
  dd.logf.sup = c(0, 2*sapply(2:(length(anchors)-1), function(i) {
      x = breaks[i+1] - anchors[i]
      (lp.breaks[i+1] - lp.anchors[i] - d.logf[i] * x)/x^2 
  }), 0)
  
  # use local polynomial approximations to build envelope
  q1 = envelope.logquad(breaks = breaks, logf = lp.anchors,
                        d.logf = d.logf,
                        dd.logf.sup = dd.logf.sup,
                        anchors = anchors)

  # draw proposal
  prop = q1$rquad(n = 1)
  
  # compute metropolis ratio (as an independence sampler)
  lR = (lp(eps = prop) - q1$dquad(x = prop, log = TRUE)) -
       (lp(eps = offset) - q1$dquad(x = offset, log = TRUE))

  # accept/reject
  if(log(runif(1)) <= lR) {
    offset = prop
    dsobs.aligned = align.obs(eps = offset)
  }
  
  
  #
  # Package results
  #

  res = list(offset = offset, dsobs.aligned = dsobs.aligned)
  
  if(debug == TRUE) {
    res$debug = list(lp = lp, q1 = q1)
  }
  
  res
}