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
#' @importFrom Matrix sparseVector Diagonal
#' 
#' @export
#' 
dsdive.impute.sample_n = function(d0, df, s0, sf, t0, tf, t.stages, rate.unif, 
                                  P.raw, P.tx, n.bins, max.tx) {
  
  # initialize output
  s.range = s0:sf
  n.stages = length(s.range)
  n = numeric(n.stages)
  
  # initialize initial state distribution
  u0 = sparseVector(x = 1, i = d0, length = n.bins)
  
  # determine window of time spent in each stage
  dt.stages = sapply(s.range, function(s) {
    min(tf, t.stages[s], na.rm = TRUE) - max(t0, t.stages[s-1])
  })
  
  # sample number of transitions in each stage
  for(i in 1:n.stages) {
    
    # extract current stage
    s = s.range[i]

    # diffuse bridging probabilities through other stages
    uf = sparseVector(x = 1, i = df, length = n.bins)
    if(i < n.stages) {
      for(j in (i+1):n.stages) {
        s2 = s.range[j]
        if(dt.stages[j] == P.raw[[s2]]$obstx.tstep) {
          uf = P.raw[[s2]]$obstx.mat %*% uf
        } else {
          uf = P.raw[[s2]]$evecs %*% 
            Diagonal(x = exp(P.raw[[s2]]$evals * dt.stages[j]), 
                     n = n.bins) %*%
            P.raw[[s2]]$evecs.inv %*% uf
        }
      }
    }
    
    #
    # compute denominator by assembling u0, uf, and within-stage transitions
    #
    
    if(dt.stages[i] == P.raw[[s]]$obstx.tstep) {
      utx = P.raw[[s]]$obstx.mat %*% uf
    } else {
      utx = P.raw[[s]]$evecs %*% 
        Diagonal(x = exp(P.raw[[s]]$evals * dt.stages[i]), n = n.bins) %*%
        P.raw[[s]]$evecs.inv %*% uf
    }
    
    C = as.numeric(u0 %*% utx)
    
    #
    # sample n[i] via inverse-transform method
    #
    
    # inverse-transform sampling variate
    uC = runif(1) * C
    
    # brute-force CDF inversion
    pC = 0
    for(n.step in 0:max.tx) {
      
      # diffuse u0 wrt. another transition
      if(n.step > 0) {
        u0 = u0 %*% P.tx[[s]]
      }
      
      # probability of reaching df after n.step transitions
      p.df = as.numeric(u0 %*% uf)
      
      # aggregate probability mass
      pC = pC + p.df * dpois(x = n.step, lambda = rate.unif * dt.stages[i])

      # stop if minimal mass is exceeded
      if(pC >= uC) {
        n[i] = n.step
        break
      }
      
      if(n.step==max.tx) {
        msg = paste('Inverse-transform sampling failed; upper bound (max.tx)',
                    'reached before inversion completed for stage', s,
                    'Excess mass:', (uC-pC)/C, sep = ' ')
        warning(msg)
      }
      
    }
  }
  
  n
}