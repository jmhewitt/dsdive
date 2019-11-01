#' Simulate dive trajectories across discrete depth bins
#'
#' The method will simulate dive trajectories from initial conditions until the 
#' trajectory is observable at \code{tf}, or a maximum number of transitions 
#' has been exceeded.  The dive simulation is bridged, so the trajectory will
#' also stop diving after returning to the surface.
#' 
#' @param depths.labels character vector that defines the depth bins
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
#'     \item{depths}{Record of which depth bins the trajectory visited}
#'     \item{durations}{Record of amount of time spent in each depth bin}
#'     \item{times}{The time at which each depth bin was entered}
#'     \item{stages}{The stage at which each depth bin was entered}
#'   }
#' 
#' @example examples/spsa.R
#' 
#' @export
#' 
#'
spsa = function(par, fn, maxit = 1e2, a = 1, c = 1, A = 1, alpha = .5, 
                gamma = 1, tol = 1e-3) {
  
  # get model dimension
  p = length(par)
  
  # optimize
  err = tol
  converged = 1
  theta = par
  ghat.last = rep(1e-6, p)
  for(k in 1:maxit) {
    
    # evaluate gain sequence
    ak = a/(A+k)^alpha
    ck = c/k^gamma
    
    # random perturbation
    delta = 2 * round(runif(p)) - 1
    
    # approximate gradient, and update parameter if finite gradient estimate
    fplus = fn(theta + ck*delta)
    if(is.finite(fplus)) {
      
      fminus = fn(theta - ck*delta)
      if(is.finite(fminus)) {
        
        # update parameter
        ghat = (fplus - fminus) / (2 * ck * delta)
        theta = theta - ak * ghat
        
        # check convergence
        err = sqrt(sum((ghat - ghat.last)^2) / sum(ghat.last^2))
        ghat.last = ghat
        if(err < tol) {
          converged = 0
          break
        }
      }
    }
  }
  
  # package results
  list(
    theta = theta,
    converged = converged
  )
}