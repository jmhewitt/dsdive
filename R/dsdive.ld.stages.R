#' Full conditional distribution for dive stages given data and parameters
#'
#' The 3 stage dive model has two breakpoints where the dive stage changes.  
#' Given data and model parameters, this function computes the conditional 
#' density of one of the breakpoints given the other.  The function is useful 
#' for constructing a Gibbs sampler that can sample over unknown stages of a 
#' dive when the dive depth bins and durations are known.  For example, this 
#' is helpful when using the \code{crawl}-based computational strategy to 
#' explore the posterior distribution for model parameters.
#' 
#' Furthermore, this full conditional distribution serves as the full 
#' conditional posterior distribution for stage breakpoints because the stage 
#' breakspoints are conditionally independent from the model parameters given 
#' the dive depth bins and durations.
#' 
#' @param breaks the two indices where a new dive stage is entered
#' @param fixed.ind  which of the stage break points in \code{breaks} should be 
#'   held fixed
#' @param beta \eqn{2 x 3} matrix in which each column contains the diving 
#'  preference and directional persistence parameters for the DIVING, SUBMERGED, 
#'  and SURFACING dive stages.
#' @param lambda length 3 vector that specifies the transition rate, 
#'   respectively in the DIVING, SUBMERGED, and SURFACING stages.
#' @param sub.tx length 2 vector that specifies the first depth bin at which 
#'   transitions to the SUBMERGED stage can occur and the probability that such 
#'   a transition occurs at the next depth transition
#' @param surf.tx parameter that specifies the probability the trajectory will 
#'   transition to the SURFACING stage at the next depth transition
#' @param depths record of depth bins the trajectory should visit
#' @param durations record of amount of time spent in each depth bin
#' @param times times at which the depth bins should be visited
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param t0.dive Time at which dive started
#'   
#' @example examples/dsdive.ld.stages.R
#' 
#' @export
#'
dsdive.ld.stages = function(breaks, fixed.ind, beta, lambda, sub.tx, surf.tx,
                            depths, durations, times, depth.bins, t0.dive) {
  
  # extract dimensions
  nt = length(times)
  
  # determine support of conditional distribution for stage breaks
  if(fixed.ind == 1) {
    if(breaks[1]+1 < nt) {
      support = (breaks[1]+1):nt
    } else {
      support = nt
    }
  } else {
    support = 1:(breaks[2]-1)
  }
  
  # compute log-mass across support
  log.wts = numeric(length(support))
  for(i in 1:length(support)) {
    # define stage vector
    if(fixed.ind == 1) {
      b = c(breaks[1], support[i])
    } else {
      b = c(support[i], breaks[2])
    }
    stage.vec = stagevec(length.out = nt, breaks = b)
    
    # compute log-weights
    log.wts[i] = dsdive.ld(depths = depths, durations = durations, 
                           times = times, stages = stage.vec, beta = beta, 
                           lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx,
                           depth.bins = depth.bins, t0.dive = t0.dive)
  }
  
  # center weights and transform to probability scale
  wts.finite = exp(scale(log.wts[is.finite(log.wts)], 
                         center = TRUE, scale = FALSE))
  
  # scale to unit mass
  wts = log.wts
  wts[is.finite(wts)] = wts.finite
  wts[is.infinite(wts)] = 0
  wts = wts / sum(wts)
  
  # return full conditional density
  data.frame(x = support, prob = wts)
}