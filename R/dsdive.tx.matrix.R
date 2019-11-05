#' Build transition matrix for dive trajectories across discrete depths
#' 
#' Dives have 3 stages, DIVING, SUBMERGED, and SURFACING.
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
#' 
#' @example examples/txmatrix.R
#' 
#' @importFrom Matrix sparseMatrix
#' 
#' @export
#' 
dsdive.tx.matrix = function(t0, depth.bins, beta, lambda, sub.tx, surf.tx) {
  
  # find maximum transition rate, to compute self-transitions
  lambda.max = max(lambda)
  
  # add a "null" depth bin to allow trajectory initialization
  num.depths = nrow(depth.bins)
  n = num.depths + 1
  
  # initialize storage for nonzero entries
  x = c()
  im = c()
  jm = c()
  
  # loop over depth bins (skip over transitions starting in "null" depth)
  for(i in 1:num.depths) {
    
    # special rules for surface bin
    if(i==1) {
    
      # all returns to the surface are absorbing
      for(s in 1:3) {
        ind = toInd(x = 2, y = 1, z = s, x.max = n, y.max = n)
        x = c(x, 1)
        im = c(im, ind)
        jm = c(jm, ind)
      }
      
      #
      # encode initial departure from surface
      #
      
      # we leave surface from stage 1, and with no last recorded depth
      p = dsdive.tx.params(t0 = t0, num.depths = num.depths, d0 = 1, 
                           d0.last = NULL, s0 = 1, beta = beta, lambda = lambda, 
                           sub.tx = sub.tx, surf.tx = surf.tx)
      
      # compute and save probability of null-transition
      self.tx = 1 - lambda[1]/lambda.max
      from.ind = toInd(x = n, y = 1, z = 1, x.max = n, y.max = n)
      if(self.tx > 0) {
        x = c(x, self.tx)
        im = c(im, from.ind)
        jm = c(jm, from.ind)
      }
      
      # transition to depth 2, stage 1
      to.ind = toInd(x = 1, y = 2, z = 1, x.max = n, y.max = n)
      prob = (1-self.tx) * (1-p$prob.stage)
      if(prob > 0) {
        x = c(x, prob)
        im = c(im, from.ind)
        jm = c(jm, to.ind)
      }
      
      # transition to depth 2, stage 2
      to.ind = toInd(x = 1, y = 2, z = 2, x.max = n, y.max = n)
      prob = (1-self.tx) * p$prob.stage
      if(prob > 0) {
        x = c(x, prob)
        im = c(im, from.ind)
        jm = c(jm, to.ind)
      }
      
    } 
    # regular depths
    else {
      
      # loop over states
      for(s in 1:3) {
        # loop over previous depths
        for(dd in c(-1,1)) {
          
          # transition parameters
          p = dsdive.tx.params(t0 = t0, num.depths = num.depths, d0 = i, 
                               d0.last = i+dd, s0 = s, beta = beta, 
                               lambda = lambda, sub.tx = sub.tx, 
                               surf.tx = surf.tx)
          
          # compute starting index
          from.ind = toInd(x = i+dd, y = i, z = s, x.max = n, y.max = n)
          
          # compute and save probability of null-transition
          self.tx = 1 - lambda[s]/lambda.max
          if(self.tx > 0) {
            x = c(x, self.tx)
            im = c(im, from.ind)
            jm = c(jm, from.ind)
          }
          
          # potential transitions to new stage
          if(s<3) {
            # transitions to new depths
            for(k in 1:length(p$labels)) {
              to.ind = toInd(x = i, y = p$labels[k], z = s+1, 
                             x.max = n, y.max = n)
              prob = (1-self.tx) * p$prob.stage * p$probs[k, s+1]
              if(prob > 0) {
                x = c(x, prob)
                im = c(im, from.ind)
                jm = c(jm, to.ind)
              }
            }
          } 
          
          # transitions to new depths within state
          for(k in 1:length(p$labels)) {
            to.ind = toInd(x = i, y = p$labels[k], z = s, 
                           x.max = n, y.max = n)
            prob = (1-self.tx) * (1-p$prob.stage) * p$probs[k, s]
            if(prob > 0) {
              x = c(x, prob)
              im = c(im, from.ind)
              jm = c(jm, to.ind)
            }
          }
        }
      }
    }
  }
  
  # build and return matrix
  nd = n^2 * 3
  m = sparseMatrix(i = im, j = jm, x = x, dims = rep(nd, 2))
  m
}