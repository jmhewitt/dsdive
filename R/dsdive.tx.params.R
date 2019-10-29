#' Compute transition parameters for dive trajectories across discrete depths
#' 
#' Dives have 3 stages, DIVING, SUBMERGED, and SURFACING.
#'   
#' @param num.depths The total number of non-surface depth bins.  Specifies the 
#'   depth bins 0, 1, ..., \code{num.depths} the model will use.  Note that the 
#'   model does not explicitly use the depth ranges the bins represent 
#'   (i.e., the bin labels) to compute transition parameters.
#' @param d0 the depth bin at which transition parameters should be computed
#' @param d0.last the previous depth bin in which the trajectory was.  If 
#'   \code{NULL}, then the autoregressive component will be skipped.
#' @param s0 the dive stage (DIVING==1, SUBMERGED==2, SURFACING==3) for which 
#'  transition parameters should be computed.  
#' @param t0 time at which transition parameters should be computed
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
#'   
#' 
#' @example examples/txparams.R
#' 
#' @export
#' 
dsdive.tx.params = function(t0, num.depths, d0, d0.last = NULL, s0, beta, 
                            lambda, sub.tx, surf.tx) {
  
  # extract probability of transitioning to next dive stage
  if(s0 == 1) { # transition to SUBMERGED stage after minimum depth reached
    if(d0 >= sub.tx[1] ) {
      prob.stage = sub.tx[2]
    } else {
      prob.stage = 0
    }
  } else if(s0 == 2) { # ASCENDING stage tx. depends on time
    prob.stage = plogis(sum(surf.tx * c(1, t0)))
  } else { # no transition out of ASCENDING stage
    prob.stage = 0
  } 
  
  # define neighboring dive bins
  if(d0 == 0) { # surface can only transition downward
    nbrs = 1
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
  res = list(rate = lambda[s0],
             prob.stage = prob.stage,
             probs = probs,
             labels = nbrs)
  
  res
}