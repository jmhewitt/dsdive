#' Simulate dive trajectories across discrete depth bins
#'
#' The method will simulate a complete dive from its beginning to end.
#' 
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param beta Directional preference model parameters.  See 
#'   \code{dsdive.tx.params} for more details.
#' @param lambda Diving rate model parameters.  See 
#'   \code{dsdive.tx.params} for more details.
#' @param t0 Initial time for the dive
#' @param steps.max maximum number of transitions to sample before stopping, 
#'   regardless of whether dive has finished.
#' @param T1 stage 1 duration
#' @param T2 stage 2 duration
#' 
#' @return A \code{dsdive} object, which is a \code{list} with the following 
#'   vectors:
#'   \describe{
#'     \item{depths}{Record of which depth bin indices the trajectory visited}
#'     \item{durations}{Record of amount of time spent in each depth bin}
#'     \item{times}{The time at which each depth bin was entered}
#'     \item{stages}{The stage at which each depth bin was entered}
#'   }
#' 
#' @example examples/dsdive.fwdsample.dive.R
#' 
#' @export
#' 
dsdive.fwdsample.dive = function(depth.bins, beta, lambda, t0, steps.max, 
                                 T1, T2) {
  
  #
  # sample stage 1 portion of dive
  #
  
  tf1 = t0 + T1
  d1 = dsdive.fwdsample.fixedstage(depth.bins = depth.bins, d0 = 1, 
                                   beta = beta, lambda = lambda, t0 = t0, 
                                   tf = tf1, steps.max = steps.max, s0 = 1)
  
  # adjust dive s.t. we stop stage 1 dynamics when stage 1 ends
  n1 = length(d1$depths)
  d1$durations[n1] = tf1 - d1$times[n1]
  
  #
  # sample stage 2 portion of dive
  #
  
  tf2 = tf1 + T2
  d2 = dsdive.fwdsample.fixedstage(depth.bins = depth.bins, 
                                   beta = beta, lambda = lambda, 
                                   d0 = d1$depths[n1],  
                                   dur0 = NULL,
                                   t0 = tf1, 
                                   tf = tf2, steps.max = steps.max, s0 = 2)
  
  
  # adjust dive s.t. we stop stage 2 dynamics when stage 2 ends
  n2 = length(d2$depths)
  d2$durations[n2] = tf2 - d2$times[n2]
  
  #
  # sample stage 3 portion of dive
  #
  
  d3 = dsdive.fwdsample.fixedstage(depth.bins = depth.bins, 
                                   beta = beta, lambda = lambda, 
                                   d0 = d2$depths[n2],  
                                   dur0 = NULL,
                                   t0 = tf2, 
                                   tf = Inf, steps.max = steps.max, s0 = 3)
  
  #
  # package dive
  #
  
  res = list(
    depths = c(d1$depths, d2$depths, d3$depths),
    times = c(d1$times, d2$times, d3$times),
    stages = c(d1$stages, d2$stages, d3$stages)
  )
  
  res$durations = c(diff(res$times), Inf)
  
  class(res) = 'dsdive'
  
  res
}