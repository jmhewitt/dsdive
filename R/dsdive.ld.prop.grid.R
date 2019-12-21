#' Likelihood for completely observed dive trajectories at grid of parameters
#'
#' Note that the function only evaluates likelihood components used for 
#' conditional densities of parameters
#' 
#' @param depths Depth bin indices in which the trajectory was observed
#' @param durations Amount of time spent in each depth bin
#' @param times Times at which each depth bin was entered
#' @param stages Record of dive stages at each depth bin
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
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param t0.dive Time at which dive started
#' @param d0.last If the depth bin that proceeded the first depth bin in 
#'   \code{depths}.  If the trajectory to be analyzed was started at the 
#'   surface, then set \code{c0.last=NULL}.
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#' 
#' @example examples/ld.R
#' 
#' @importFrom stats dexp dbinom
#' 
#' @param pi1 vector of grid points for pi1 downward probability
#' @param pi2 vector of grid points for pi2 downward probability
#' 
#' @export
#' 
dsdive.ld.prop.grid = function(depths, durations, stages, pi1, pi2, lambda1, 
                               lambda2, lambda3, depth.bins) {
  
  pi1.len = length(pi1)
  pi2.len = length(pi2)
  lambda1.len = length(lambda1)
  lambda2.len = length(lambda2)
  lambda3.len = length(lambda3)
  
  #
  # extract dive components
  #
  
  # stage filters
  s1 = stages==1
  s2 = stages==2
  s3 = stages==3
  
  # durations split by stages
  dur.stage.1 = durations[s1]
  dur.stage.2 = durations[s2]
  dur.stage.3 = durations[s3]
  # depths split by stages
  dep.stage.1 = depths[s1]
  dep.stage.2 = depths[s2]
  dep.stage.3 = depths[s3]
  # depth bin widths split by stages
  width.stage.1 = 2 * depth.bins[dep.stage.1, 2]
  width.stage.2 = 2 * depth.bins[dep.stage.2, 2]
  width.stage.3 = 2 * depth.bins[dep.stage.3, 2]
  
  # total transitions in each stage
  n.1 = sum(s1)
  n.3 = sum(s3)
  # downward transitions in each stage
  inds.1 = which(s1)
  inds.3 = which(s3)
  down.1 = sum(diff(depths[min(inds.1):(max(inds.1)+1)]) == 1)
  if(n.3 > 0) {
    down.3 = sum(diff(depths[min(inds.3):(max(inds.3))]) == 1)
  } else {
    down.3 = 0
  }
  
  
  
  #
  # evaluate density components along grid dimensions
  #
  
  # log-densities of durations
  ld.dur.1 = sapply(lambda1, function(lambda) {
    sum(dexp(x = dur.stage.1, rate = lambda / width.stage.1, log = TRUE))
  })
  ld.dur.2 = sapply(lambda2, function(lambda) {
    sum(dexp(x = dur.stage.2, rate = lambda / width.stage.2, log = TRUE),
        na.rm = TRUE)
  })
  if(n.3 > 0) {
    ld.dur.3 = sapply(lambda3, function(lambda) {
      sum(dexp(x = dur.stage.3, rate = lambda / width.stage.3, log = TRUE), 
          na.rm = TRUE)
    })
  } else {
    ld.dur.3 = rep(0, lambda3.len)
  }
  
  
  # log-densities of depth bin transitions
  ld.tx.1 = sapply(pi1, function(p) {
    dbinom(x = down.1, size = n.1, prob = p, log = TRUE)
  })
  if(n.3 > 0) {
    ld.tx.2 = sapply(pi2, function(p) {
      dbinom(x = down.3, size = n.3, prob = p, log = TRUE)
    })
  } else {
    ld.tx.2 = rep(0, lambda3.len)
  }
  
  
  #
  # assemble log-density across grid
  #
  
  ld = array(data = 0, 
             dim = c(pi1.len, pi2.len, lambda1.len, lambda2.len, lambda3.len))
  
  for(i in 1:pi1.len) {
    li = ld.tx.1[i]
    for(j in 1:pi2.len) {
      lj = ld.tx.2[j]
      for(k in 1:lambda1.len) {
        lk = ld.dur.1[k]
        for(l in 1:lambda2.len) {
          ll = ld.dur.2[l]
          for(m in 1:lambda3.len) {
            ld[i,j,k,l,m] = li + lj + lk + ll + ld.dur.3[m]
  }}}}}
  
  as.numeric(ld)
}