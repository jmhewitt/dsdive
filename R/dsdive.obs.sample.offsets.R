#' Sample stage transition times from full conditional posterior
#'
#' We assume the prior distributions for the stage transition times are Gamma
#' distributions with shape/rate parameterization.
#'
#' @example examples/dsdive.obs.sample.offsets.R
#'
#'
dsdive.obs.sample.offsets = function(dsobs.aligned, dsobs.unaligned, offset,
                                     offset.tf, t.stages, P.raw, depth.bins, 
                                     tstep, max.width, prior.params,
                                     sample.start, debug = FALSE) {

  tstep2 = 2*tstep

  #
  # components for constructing proposal densities for full conditionals
  #

  lp = function(eps) { sapply(eps, function(eps) {
    # log-posterior for data given a dive offset
    #
    # Parameters:
    #  eps - offset
    
    # set offsets depending on whether start or end of dive is being sampled
    if(sample.start) {
      eps.t0 = eps
      eps.tf = offset.tf
    } else {
      eps.t0 = offset
      eps.tf = eps
    }
    
    # construct dsobs object with realigned observation times
    aligned = dsdive.align.obs(depths = dsobs.unaligned$depths,
                               times = dsobs.unaligned$times,
                               t.stages = t.stages, offset = eps.t0, 
                               offset.tf = eps.tf)

    # log density associated with dive offset
    ld = dsdive.obsld(dsobs.list = aligned, t.stages.list = t.stages,
                      P.raw = P.raw, s0 = 1, sf = 3)

    # return log-posterior (ignoring jacobian)
    ld + dbeta(x = (eps + tstep)/tstep2,
               shape1 = prior.params[1], shape2 = prior.params[2],
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

  # determine slope and curvature for each polynomial segment of envelope
  envelope = sapply(1:length(anchors), function(i) {
    # extract time-coordinates for interval
    L.x = breaks[i]
    M.x = anchors[i]
    U.x = breaks[i+1]
    # extract log-posterior for interval
    L = lp.breaks[i]
    M = lp.anchors[i]
    U = lp.breaks[i+1]
    # the in/finiteness of the points determines processing
    L.undef = !is.finite(L)
    M.undef = !is.finite(M)
    U.undef = !is.finite(U)

    # assume lp = -Inf throughout entire interval
    if(L.undef & M.undef & U.undef) {
      c(0, 0, -Inf)
    }
    # assume lp = -Inf at start and midpoint, but ends finite
    else if(L.undef & M.undef) {
      c(0, 0, U)
    }
    # assume lp is only -Inf at start point
    else if(L.undef) {
      c(0, (U-M)/(U.x-M.x), M)
    }
    # assume lp is -Inf at end and midpoint, but starts finite
    else if(U.undef & M.undef) {
      c(0, 0, L)
    }
    # assume lp is only -Inf at end point
    else if(U.undef) {
      c(0, (L-M)/(L.x-M.x), M)
    }
    # assume lp is finite in entire interval
    else {
      b = (U-L)/(U.x-L.x)
      x = U.x - M.x
      a = (U - M - b * x)/x^2
      c(2*a, b, M)
    }
  })

  # extract slopes and curvatures from envelope
  d.logf = envelope[2,]
  dd.logf.sup = envelope[1,]

  # zero-out small curvatures, which are numerically unstable in logquad fn.
  dd.logf.sup[abs(dd.logf.sup)<1e-4] = 0

  # use local polynomial approximations to build envelope
  q1 = envelope.logquad(breaks = breaks, logf = envelope[3,],
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
    
    # set offsets depending on whether start or end of dive is being sampled
    if(sample.start) {
      eps.t0 = eps
      eps.tf = offset.tf
    } else {
      eps.t0 = offset
      eps.tf = eps
    }
    
    # construct dsobs object with realigned observation times
    dsobs.aligned = dsdive.align.obs(depths = dsobs.unaligned$depths,
                                     times = dsobs.unaligned$times,
                                     t.stages = t.stages, offset = eps.t0,
                                     offset.tf = eps.tf)
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
