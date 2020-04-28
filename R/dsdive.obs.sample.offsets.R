#' Sample t0 and tf offsets from full conditional posteriors
#'
#' Sampler uses a piecewise quadratic polynomial approximation to the log of the 
#' full conditional posterior.  Note that this sampler will only sample one 
#' offset at a time.
#'
#' @param dsobs.aligned \code{dsobs} object that has been adjusted to account 
#'   for the t0 and tf offsets.  Note that this should be based on the current 
#'   value of the offsets, and that this function will update these parameters.
#' @param dsobs.unaligned \code{dsobs} object with the raw observations of a 
#'   dive object
#' @param offset value of current t0 offset
#' @param offset.tf value of current tf offset
#' @param t.stages stage transition times for the dive
#' @param P.raw list of continuous time probability transition matrices, and 
#'  components.
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param tstep Time between observations in \code{dsobs.unaligned}
#' @param max.width The t0 and tf offsets are updated with a piecewise 
#'   proposal distribution.  This parameter controls the maximum width 
#'   of the intervals for the proposal distribution.  This is a tuning parameter 
#'   that controls the numerical stability of the proposal distribution, which 
#'   is sampled via inverse CDF techniques.
#' @param prior.params \code{shape1} and \code{shape2} parameters for the 
#'   scaled and shifted Beta prior distribution for the offset being sampled.
#' @param sample.start \code{TRUE} to sample t0 offset; 
#'   \code{FALSE} to sample tf offset.
#' @param debug \code{TRUE} to return debugging objects, such as the proposal 
#'   density and log posterior
#' @param lpapprox (Optional) Object containing the piecewise-quadratic
#'   approximation to the log-posterior used to propose model parameters.  
#'   If \code{NULL}, then an approximation will be computed.
#' @param output.lpapprox \code{TRUE} to return the approximation used 
#'   to propose model parameters.
#'   
#' # @example examples/dsdive.obs.sample.offsets.R
#' 
#' @importFrom stats dbeta runif
#' 
#' @export
#'
dsdive.obs.sample.offsets = function(dsobs.aligned, dsobs.unaligned, offset,
                                     offset.tf, t.stages, P.raw, depth.bins, 
                                     tstep, max.width, prior.params,
                                     sample.start, debug = FALSE, 
                                     lpapprox = NULL, output.lpapprox = FALSE) {

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
  # build or extract proposal distribution
  #
  
  if(!is.null(lpapprox)) {
    
    q1 = lpapprox
    
  } else {
    
    # determine intervals and midpoints for log-quadratic sampling envelope
    breaks = refine.partition(
      breaks = c(-tstep, 0, tstep),
      max.width = max.width)
    anchors = breaks[1:(length(breaks)-1)] + diff(breaks)/2
    
    # eval log-posterior at all points; account for Beta prior boundary issues
    tol.boundary = .Machine$double.eps * tstep2
    end.inds = c(1,length(breaks))
    lp.breaks = lp(eps = c(breaks[1] + tol.boundary,
                           breaks[-end.inds],
                           breaks[length(breaks)] - tol.boundary))
    lp.anchors = lp(eps = anchors)
    
    # determine slope and curvature for each polynomial segment of envelope
    envelope = envelope.approx(breaks = breaks, anchors = anchors, 
                               lp.breaks =lp.breaks, lp.anchors = lp.anchors)
    
    # zero-out small curvatures, which are numerically unstable in logquad fn.
    envelope$dd.logf.sup[abs(envelope$dd.logf.sup)<1e-4] = 0
    
    # use local polynomial approximations to build envelope
    q1 = envelope.logquad(breaks = breaks, logf = envelope$logf,
                          d.logf = envelope$d.logf,
                          dd.logf.sup = envelope$dd.logf.sup,
                          anchors = anchors)
  }
  

  #
  # sample dive offset
  #
  
  # draw proposal
  prop = q1$rquad(n = 1)

  # compute metropolis ratio (as an independence sampler)
  lR = (lp(eps = prop) - q1$dquad(x = prop, log = TRUE)) -
       (lp(eps = offset) - q1$dquad(x = offset, log = TRUE))

  # accept/reject
  accept = log(runif(1)) <= lR
  if(accept) {
    
    if(sample.start) {
      offset = prop
    } else {
      offset.tf = prop
    }
    
    # construct dsobs object with realigned observation times
    dsobs.aligned = dsdive.align.obs(depths = dsobs.unaligned$depths,
                                     times = dsobs.unaligned$times,
                                     t.stages = t.stages, offset = offset,
                                     offset.tf = offset.tf)
  }


  #
  # Package results
  #

  res = list(offset = ifelse(sample.start, offset, offset.tf), 
             dsobs.aligned = dsobs.aligned, accepted = accept)

  if(debug == TRUE) {
    res$debug = list(lp = lp, q1 = q1)
  }
  
  if(output.lpapprox) { 
    res$q1 = q1
  }

  res
}
