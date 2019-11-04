#' Simulate dive trajectories across discrete depth bins
#'
#' The method will simulate dive trajectories from initial conditions until the 
#' trajectory is observable at \code{tf}, or a maximum number of transitions 
#' has been exceeded.  The dive simulation is bridged, so the trajectory will
#' also stop diving after returning to the surface.
#' 
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param d0 the depth bin at which transition parameters should be computed
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
#' @param t0 time at which transition parameters should be computed
#' @param tf time at which sampling should end after
#' @param steps.max maximum number of transitions to sample before stopping, 
#'   regardless of whether \code{tf} is reached.
#' @param dur0 time spent at location \code{d0}.  If \code{NULL}, then a 
#'   duration in state \code{d0} will be sampled, otherwise a new state will 
#'   be sampled first, then sampling will continue from the new state at time 
#'   \code{t0 + dur0}.
#' @param s0 dive stage in which forward simulation begins
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
#' @example examples/txparams.R
#' 
#' @export
#' 
#'
dsdive.fwdsample = function(depths.labels, d0, beta, lambda, sub.tx, surf.tx, 
                            t0, tf, steps.max, dur0 = NULL, nsteps = NULL, s0) {
  
  # extract information about discretized depth domain
  num.depths = nrow(depths.labels) - 1
  
  # initialize output components
  depths = d0
  durations = dur0
  times = t0
  stages = s0
  
  # validate time window
  if(t0 > tf) {
    stop('Dive sampling only goes forward in time, but initial time t0 is past 
         target time tf.')
  }
  
  # sample stages and transitions if trajectory is not observable at tf
  if(t0 + ifelse(is.null(dur0), 0, dur0) < tf) {
    
    # initialize state tracking
    current = list(
      depth = d0,
      duration = dur0,
      t = t0,
      depth.last = NULL,
      stage = stages
    )
    
    # configure so that we sample an exact number of transitions
    if(!is.null(nsteps)) {
      tf = Inf
      steps.max = nsteps + 1 
    }
    
    # forward sample up to maximum number of steps
    for(i in 1:steps.max) {
      
      # get sampling parameters for this location
      params.tx = dsdive.tx.params(t0 = current$t, num.depths = num.depths,
                                   d0 = current$depth, s0 = current$stage,
                                   d0.last = current$depth.last,  beta = beta, 
                                   lambda = lambda, sub.tx = sub.tx, 
                                   surf.tx = surf.tx)
      
      # if necessary, sample and save duration for current state
      if(is.null(current$duration)) {
        dur = rexp(n = 1, rate = params.tx$rate) 
        durations = c(durations, dur)
      } else {
        dur = current$duration
      }
      
      # sample stage
      stage = current$stage + rbinom(1, 1, params.tx$prob.stage)
      # sample new depth
      if(length(params.tx$labels)==1) {
        d = params.tx$labels
      } else {
        d = sample(x = params.tx$labels, size = 1, 
                   prob = params.tx$probs[, stage])
      }
      
      # update state to prep for next transition
      current = list(
        depth = d,
        duration = NULL, # duration in new state is not yet known
        t = current$t + dur,
        depth.last = current$depth,
        stage = stage
      )
      
      # append new state and its transition time to output
      depths = c(depths, d)
      times = c(times, current$t)
      stages = c(stages, stage)
      
      # stop sampling once trajectory is observable at tf
      if(current$t >= tf)
        break
      # or stop sampling once trajectory returns to surface
      else if(current$depth == 0)
        break
    }
    
    # only retain samples so that the trajectory is observable at tf
    times.observable = which(times>=tf)
    if(any(times.observable)) {
      inds = 1:(min(times.observable) - 1)
    } else {
      if(current$depth == 0) {
        inds = 1:length(depths)
      } else if(!is.null(nsteps)) {
        inds = 1:(nsteps+1)
      } else {
        inds = NULL
      }
    }
    
    depths = depths[inds]
    durations = durations[inds]
    times = times[inds]
    stages = stages[inds]
  }
  
  
  #
  # package results
  #
  
  res = list(
    depths = depths,
    durations = durations,
    times = times,
    stages = stages
  )
  
  class(res) = 'dsdive'
  
  res
}