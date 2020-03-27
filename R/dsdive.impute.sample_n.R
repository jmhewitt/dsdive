#' Sample number of transitions between observations for a CTMC
#' 
#' Uses properties of homogeneous Continuous time Markov Chains (CTMCs) to 
#' sample the number of transitions between known endpoints and times.
#'   
#' @param d0 index of depth bin in which the CTMC segment starts
#' @param df index of depth bin in which the CTMC segment ends
#' @param s0 dive stage in which the CTMC begins
#' @param sf dive estage in which the CTMC ends
#' @param t0 time at which \code{d0} is observed
#' @param tf time at which \code{df} is observed
#' @param t.stages stage transition times for the dive
#' @param rate.unif uniformization rate, for standardizing transition
#'   rates between states
#' @param P.raw list of continuous time probability transition matrices, and 
#'  components.
#' @param P.tx list of discrete time probability transition matrices
#' @param n.bins number of rows in the \code{depths} matrix
#' @param max.tx maximum number of transitions between observations that will 
#'   be allowed during imputation
#' 
#' @importFrom Matrix sparseVector Diagonal
#' @importFrom stats runif dpois
#' @importFrom expm expAtv
#' 
# @example examples/dsdive.impute.sample_n.R
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
          uf = expAtv(A = as.matrix(P.raw[[s2]]$A), 
                      t = dt.stages[j],
                      v = as.numeric(uf))[[1]]
        }
      }
    }
    
    #
    # compute denominator by assembling u0, uf, and within-stage transitions
    #
    
    if(dt.stages[i] == P.raw[[s]]$obstx.tstep) {
      utx = P.raw[[s]]$obstx.mat %*% uf
    } else {
      utx = expAtv(A = as.matrix(P.raw[[s]]$A), 
                   t = dt.stages[i],
                   v = as.numeric(uf))[[1]]
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