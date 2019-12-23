#' Compute stage transition parameters for dive trajectories across discrete depths
#'
#' Computes the probability that a dive will transition between any of the 3 
#' stages, PRIMARY DESCENT (PD), INTERMEDIATE BEHAVIORS (IB), and 
#' PRIMARY ASCENT (PA).
#'   
#' @param t0 vector of times at which transition probs. should be computed
#' @param sub.tx length 2 vector that specifies the first depth bin at which 
#'   transitions to the IB stage can occur and the probability that such 
#'   a transition occurs at the next depth transition
#' @param surf.tx parameter that specifies the probability the trajectory will 
#'   transition to the PA stage at the next depth transition
#' @param t0.dive Time at which dive started
#' @param t.scale Nominally, \code{t0} and \code{t0.dive} should represent times 
#'   in seconds relative to a starting epoch (e.g., 0 for simulations, or 
#'   1/1/1970 for UTC date), however it may make sense to define the stage 
#'   transition model parameters with respect to a different time scale 
#'   (e.g., minutes, hours).  The argument \code{t.scale} will rescale \code{t0}
#'   and \code{t0.dive} to accomplish this.
#' @param t.stage2 time at which second stage was entered.  If stage 2 has not 
#'   yet been entered, then set \code{t.stage2 = NA}, and 
#'   \code{dsdive.tx.stage} will implicitly set \code{t.stage2 = t0}.
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#' 
#' @example examples/dsdive.tx.stage.R
#' 
#' @importFrom stats pexp
#' 
#' @export
#' 
dsdive.tx.stage = function(t0, s0, sub.tx, surf.tx, t0.dive, t.stage2,
                           t.scale = 60, rates, model) {
  
  if(is.na(t.stage2)) {
    t.stage2 = t0
  }
  
  #
  # compute probabilities
  #
  
  if(model == 'conditional') {
    
    # probability of transition to stage 2
    if(s0 == 1) {
      sub.lu = c(sub.tx[1] - sub.tx[2], sub.tx[1] + sub.tx[2])
      tgt.time = sub.lu * t.scale + t0.dive
      
      interval.start = ifelse(t0 < tgt.time[1], tgt.time[1], t0)
      if(interval.start < tgt.time[2]) {
        prob.inrange = diff(pexp(q = c(interval.start, tgt.time[2]) - t0, 
                                 rate = rates))
        interval = max(tgt.time[2] - interval.start, 0)
        scale.prob = 1 / (rates * interval)
        res = min(prob.inrange * scale.prob, 1)
      } else {
        res = 1
      }
    } 
    # probability of transition to stage 3
    else if(s0 == 2) {
      surf.lu = c(surf.tx[1] - surf.tx[2], surf.tx[1] + surf.tx[2])
      tgt.time = surf.lu * t.scale + t.stage2
      
      interval.start = ifelse(t0 < tgt.time[1], tgt.time[1], t0)
      if(interval.start < tgt.time[2]) {
        prob.inrange = diff(pexp(q = c(interval.start, tgt.time[2]) - t0, 
                                 rate = rates))
        interval = max(tgt.time[2] - interval.start, 0)
        scale.prob = 1 / (rates * interval)
        res = min(prob.inrange * scale.prob, 1)
      } else {
        res = 1
      }
    }
    # probability of transition out of stage 3 is 0
    else if(s0 == 3) {
      res = 0
    }
    
  } else if(model == 'logit') {
    
    # y = qlogis(p = c(.01,.99))
    y = c(-4.59512,  4.59512)
    
    if(s0 < 3) {
      
      if(s0 == 1) {
        sub.lu = c(sub.tx[1] - sub.tx[2], sub.tx[1] + sub.tx[2])
        tgt.time = sub.lu * t.scale + t0.dive
      } else if(s0 == 2) {
        surf.lu = c(surf.tx[1] - surf.tx[2], surf.tx[1] + surf.tx[2])
        tgt.time = surf.lu * t.scale + t.stage2
      }
      
      beta2 = diff(y)/diff(tgt.time)
      beta1 = y[1] - beta2 * tgt.time[1]
      res = plogis(beta1 + beta2 * t0)
      
    } else if(s0 == 3) {
      res = 0
    }
  
  }
   
  res
}