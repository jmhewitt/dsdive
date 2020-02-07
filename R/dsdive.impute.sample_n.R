#' Compute transition matrix for Markov chain representation of dive model
#' 
#' Given model parameters and transition times, the 3 stage dive model can be 
#' represented as a Markov chain over a state space where each state records
#' the last depth bin, current depth bin, and current dive stage.  This 
#' function computes the complete transition matrix for the state space, given 
#' model parameters and time at which the transition will take place.
#'   
#' @param t0 time at which transition parameters should be computed
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
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
#' @param inflation.factor.lambda In order to facilitate bridged transitions, 
#'   the transition rate of the overall process must be inflated to allow the 
#'   possibility of self-transitions.  Self-transitions allow bridged paths to 
#'   dynamically modify the total number of transitions between observed values
#'   so that a valid path between observations is always possible.  The 
#'   \code{inflation.factor.lambda} parameter implicitly controls the number of 
#'   self-transitions that will occur.  Larger values will create more 
#'   self-transitions.
#' @param min.depth As a computational efficiency option, only compute 
#'   transition parameters for depth bins at and above \code{min.depth}.
#' @param max.depth As a computational efficiency option, only compute 
#'   transition parameters for depth bins at and below \code{max.depth}.
#' @param t0.dive Time at which dive started
#' @param lambda.max Arrival rate for the parent Poisson process that will
#'   be thinned.  \code{lambda.max} will be scaled by 
#' @param t.stage2 time at which second stage was entered
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#'   
#' @example examples/dsdive.impute.sample_n.R
#' 
#' @param n0 number of transitions during initial stage; if NULL, then this will 
#'  be the number returned.  Otherwise, nf, the number of transitions during 
#'  final stage will be returned instead.
#'
#' @param ff.s0 list of n-step forward distributions when a DTMC is started at 
#'  node d0.
#' @param P.raw list of continuous time probability transition matrices, and 
#'  components.
#' @param P.tx list of discrete time probability transition matrices
#' 
#' @export
#' 
dsdive.impute.sample_n = function(
  n0=NULL, d0, df, s0, sf, t0, tf, t.stages, rate.unif, P.raw, P.tx, 
  ff.s0 = NULL, ff.sf = NULL, n.bins, max.tx, ff.out = FALSE) {
  
  #
  # initialize n-step forward distributions if necessary
  #
  
  if(is.null(ff.s0)) {
    ff.s0 = vector('list', length = 0)
    a = numeric(length = n.bins)
    a[d0] = 1
    ff.s0[[1]] = a
  }
  
  if(is.null(ff.sf)) {
    ff.sf = vector('list', length = 0)
    a = numeric(length = n.bins)
    a[df] = 1
    ff.sf[[1]] = a
  }
  
  
  #
  # compute normalizing constant, sampling window, and transition constants
  #
  
  if(s0 == sf) {
    # sampling window
    twin = tf - t0
    # within-stage probability of observing the given transition
    C = P.raw[[s0]]$obstx.mat[d0,df]
  } else {
    # depth bin distribution from stage transition to next observation
    uf = P.raw[[sf]]$evecs %*% 
      (exp(P.raw[[sf]]$evals * (tf-t.stages[sf-1])) * 
         P.raw[[sf]]$evecs.inv[,df]) 
    # depth bin distribution from d0 to stage transition
    if(is.null(n0)) {
      # sampling window
      twin = t.stages[sf-1] - t0
      # number of transitions is unknown
      u0 = t((P.raw[[s0]]$evecs[d0,] *  exp(P.raw[[s0]]$evals * twin)) %*% 
             P.raw[[s0]]$evecs.inv)
    } else {
      # sampling window
      twin = tf - t.stages[sf-1]
      # update n-step forward distributions if necessary
      if(length(ff.s0) <= n0) {
        # number of additional transitions needed
        steps = n0 - length(ff.s0) + 1
        # all transitions are to unknown locations
        L = matrix(1/n.bins, nrow = n.bins, ncol = steps)
        # transition matrices for each step
        B = lapply(1:steps, function(i) P.tx[[s0]])
        # append additional forward distributions
        ff.s0 = c(ff.s0, ff(B = B, L = L, a0 = ff.s0[[length(ff.s0)]])[-1])
      }
      # number of transitions is known
      u0 = ff.s0[[length(ff.s0)]]
    }
    # between-stage probability of observing the given transition
    C = as.numeric(t(u0) %*% uf)
  }
  
  
  #
  # sample n via inverse-transform method
  #
  
  # inverse-transform sampling variate
  uC = runif(1) * C
  
  # brute-force CDF inversion
  pC = 0
  for(n in 0:max.tx) {
    
    # update n-step forward distributions if necessary
    if(is.null(n0)) {
      if(length(ff.s0) <= n) {
        ff.s0[[n+1]] = t(P.tx[[s0]]) %*% ff.s0[[n]]
      }
    } else {
      if(length(ff.sf) <= n) {
        ff.sf[[n+1]] = t(P.tx[[sf]]) %*% ff.sf[[n]]
      }
    }
    
    # probability of reaching df after n transitions
    if(s0==sf) {
      # within-stage transition
      p.df = ff.s0[[n+1]][df]
    } else {
      if(is.null(n0)) {
        # between stage transition and number of s0 transitions is unknown
        p.df = as.numeric(t(ff.s0[[n+1]]) %*% uf)
      } else {
        # between stage transition and number of s0 transitions is known
        p.df = as.numeric(t(u0) %*% ff.sf[[n+1]])
      }
    }
    
    # aggregate probability mass
    pC = pC + p.df * dpois(x = n, lambda = rate.unif * twin)
    
    # stop if minimal mass is exceeded
    if(pC >= uC) {
      break
    }
  }
  
  if(n==max.tx) {
    msg = paste('Inverse-transform sampling failed; upper bound (max.tx)',
                'reached before inversion completed.  Excess mass:',
                (uC-pC)/C, sep = ' ')
    warning(msg)
  }
  
  if(ff.out) {
    list(n = n, ff.s0 = ff.s0, ff.sf = ff.sf)
  } else {
    n
  }
}