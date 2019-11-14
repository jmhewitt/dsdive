#' Compute transition parameters for dive trajectories across discrete depths
#'
#' Computes elements of a transition matrix when given model parameters and 
#' time/space locations for a dive model that has 3 stages, PRIMARY DESCENT 
#' (PD), INTERMEDIATE BEHAVIORS (IB), and PRIMARY ASCENT (PA).
#'   
#' @param t0 time at which transition parameters should be computed
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param d0 the depth bin at which transition parameters should be computed
#' @param d0.last the previous depth bin in which the trajectory was.  If 
#'   \code{NULL}, then the autoregressive component will be skipped.
#' @param s0 the dive stage (PRIMARY DESCENT==1, INTERMEDIATE BEHAVIORS==2, 
#'   PRIMARY ASCENT==3) for which transition parameters should be computed.  
#' @param beta \eqn{2 x 3} matrix in which each column contains the diving 
#'  preference and directional persistence parameters for the PRIMARY DESCENT 
#'  (PD), INTERMEDIATE BEHAVIORS (IB), and PRIMARY ASCENT (PA) dive stages.
#' @param lambda length 3 vector that specifies the transition rate, 
#'   respectively in the PD, IB, and PA stages.
#' @param sub.tx length 2 vector that specifies the first depth bin at which 
#'   transitions to the IB stage can occur and the probability that such 
#'   a transition occurs at the next depth transition
#' @param surf.tx parameter that specifies the probability the trajectory will 
#'   transition to the PA stage at the next depth transition
#' @param t0.dive Time at which dive started
#' 
#' @example examples/txparams.R
#' 
#' @export
#' 
dsdive.tx.params = function(t0, depth.bins, d0, d0.last = NULL, s0, beta, 
                            lambda, sub.tx, surf.tx, t0.dive) {
  
  num.depths = nrow(depth.bins)
  
  # extract probability of transitioning to next dive stage
  prob.stage = dsdive.tx.stage(t0 = t0, d0 = d0, sub.tx = sub.tx, 
                               surf.tx = surf.tx, t0.dive = t0.dive)[s0,]
  
  # define neighboring dive bins
  if(d0 == 1) { # surface can only transition downward
    nbrs = 2
  } else if(d0 == num.depths) { # max depth can only transition upward
    nbrs = num.depths - 1
  } else { # all other depths can go up or down one bin
    nbrs = d0 + c(-1, 1)
  }
  
  # define transition probabilities between dive bins for each dive stage
  if(length(nbrs) == 1) {
    probs = rep(1,3)
  } else {
    # diving preference component
    logit.probs = matrix(c(-beta[1,], beta[1,]), nrow = 2, byrow = TRUE)
      
    # add autoregressive component
    if(any(beta[2,] != 0)) {
      if(!is.null(d0.last)) {
        # define directions to neighbors
        dir.nbrs = c(-1,1)
        # compute direction to last node
        dir.last = dir.nbrs[which(d0.last==nbrs)]
        # scale direction by relative bin sizes
        dir.last = dir.last * depth.bins[d0.last, 2] / depth.bins[d0, 2]
        # add autoregressive component
        dirs.ar = dir.nbrs * dir.last
        for(i in 1:3) {
          logit.probs[,i] = logit.probs[,i] + beta[2,i] * dirs.ar
        }
      }
    }
    # convert probabilities
    probs = plogis(logit.probs)
  }
  
  
  # package results
  res = list(rate = lambda[s0] / (2 * depth.bins[d0, 2]),
             prob.stage = prob.stage,
             probs = matrix(probs, ncol = 3),
             labels = nbrs)
  
  res
}