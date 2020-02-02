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
#' @param t.stage2 time at which second stage was entered
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#'   
#'
#' @importFrom Matrix sparseVector
#' 
#' @example examples/dsdive.obsld.R
#' 
#' @export
#'
dsdive.obsld = function(dsobs.list, t.stages.list, P.raw) {
  
  # number of depth bins
  n.bins = nrow(P.raw[[1]]$obstx.mat)
  
  # compute log-density for family of observations
  sum(sapply(1:length(dsobs.list), function(dive.id) {
    
    # extract observed dive components
    d = dsobs.list[[dive.id]]
    depths = d$depths
    times = d$times
    n = length(d$depths)
    
    # compute stage at each observation
    t.stages = t.stages.list[[dive.id]]
    stages = findInterval(times, t.stages) + 1
    
    ld = 0
    for(i in 1:(n-1)) {
      # extract start/end depth bins
      d0 = depths[i]
      df = depths[i+1]
      # extract start/end stages
      s0 = stages[i]
      sf = stages[i+1]
      # add likelihood contribution of observation
      if(s0==sf) { 
        # add contribution for within-stage transition
        ld = ld + log(P.raw[[s0]]$obstx.mat[d0,df])
      } else {
        
        #
        # build contribution for between-stage transitions as a series of 
        # matrix-vector products
        #
        
        # initial depth bin
        u0 = sparseVector(x = 1, i = d0, length = n.bins)
        # depth bin distribution at time of stage transition
        u0 = ((t(u0) %*% P.raw[[s0]]$evecs) * 
              exp(P.raw[[s0]]$evals * (t.stages[s0] - times[i]))) %*% 
             P.raw[[s0]]$evecs.inv
        
        # final depth bin
        uf = sparseVector(x = 1, i = df, length = n.bins)
        # depth bin distribution at end of stage transition
        uf = ((t(uf) %*% t(P.raw[[sf]]$evecs.inv)) * 
              exp(P.raw[[sf]]$evals * (times[i+1] - t.stages[s0]))) %*% 
             t(P.raw[[sf]]$evecs)
        
        # add contribution for between-stage transition
        ld = ld + log(sum(u0 * uf))
      }
    }
    
    ld
  }))
  
}