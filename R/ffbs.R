#' Use bridged sampling to impute a complete dive trajectory consistent with observations
#'
#' The sampling method is designed to sample many trajectories simultaneously, 
#' so has an extra level of approximation in the proposal distributions.
#' 
#' @param M the number of trajectories to sample
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param depths record of depth bins the trajectory should visit
#' @param times times at which the depth bins should be visited
#' @param s0 dive stage at which the trajectory should be started from
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
#' @param verbose If \code{TRUE}, then the sampler's progress will be printed 
#'   during sampling.
#' @param precompute.bridges If \code{TRUE}, then the bridged transition 
#'   matrices will be precomputed.  Enabling this option will increase the 
#'   memory overhead of the method, but will reduce its runtime.
#' @param t0.dive Time at which dive started
#' @param resample Resample particles at each step if \code{TRUE}.
#' @param trajectory.conditional If not \code{NULL}, then 
#'   \code{trajectory.conditional} must be the dive information for a 
#'   completely observed \code{dsdive} object.  The first entries in the 
#'   initialization vectors \code{d0}, \code{d0.last}, \code{s0} must be 
#'   associated with the trajectory observed in \code{trajectory.conditional}.
#'   Providing a non \code{NULL} value for \code{trajectory.conditional} will 
#'   cause \code{dsdive.fastbridge} to simulate one fewer trajectories, as 
#'   the value of \code{trajectory.conditional} will be returned.  The 
#'   importance of this function argument is that \code{dsdive.fastbridge}
#'   will evaluate the proposal density for \code{trajectory.conditional}.
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#' @param ld.compute \code{TRUE} to compute likelihood values as well.  This 
#'   is required if resampling or conditional trajectory imputation is used.
#' @param t.stages the times at which stage transitions occur
#' 
#' @example examples/ffbs.R
#' 
#' @export
#'
ffbs = function(Bt, O, W, lambda, shift.scale = 0) {
  
  # extract initial conditions
  s.init = O[1]
  
  # get size of state space
  N = nrow(Bt[[1]])
  
  # get number of observations
  T = length(O)
  
  # build likelihood matrix; NA entries are non-informative
  L = matrix(0, nrow = N, ncol = T)
  for(i in 1:T){
    if(!is.na(O[i])) {
      L[O[i],i] = 1
    } else {
      L[,i] = 1/N
    }
  }
  
  # initialize forward-filtering vectors
  a = vector('list', T+1)
  
  # encode initial state distribution
  a[[1]] = numeric(N)
  a[[1]][s.init] = 1
  
  # forward filter
  for(i in 2:length(a)) {
    # at = numeric(N)
    # for(j in 1:N) {
    #   at[j] = sum(a[[i-1]] * L[,i-1] * Bt[[i-1]][,j])
    # }
    # a[[i]] = at
    a[[i]] = t(Bt[[i-1]]) %*% (a[[i-1]] * L[,i-1])
    # rescale a
    a[[i]] = a[[i]] / sum(a[[i]])
    if(any(is.na(as.numeric(a[[i]])))) {
      browser()
    }
  }
  
  # browser()
  
  # backward sample
  s = numeric(T)
  p = as.numeric(L[,T] * a[[T+1]])

  p.max = which.max(p)
  p[p.max] = max(shift.scale, p[p.max])
  p[-p.max] = p[-p.max]/sum(p[-p.max]+.Machine$double.eps) * (1-p[p.max])

  
  if(any(is.na(p))) {
    browser()
  }
  
  
  
  s[T] = sample(x = 1:N, size = 1, prob = p)
  # s[T] = which.max(p)
  for(t in (T-1):1) {
    # s[t] = sample(x = 1:N, size = 1, prob = a[[t]] * Bt[[t]][,s[t+1]] * L[,t])
    p = as.numeric(a[[t]] * Bt[[t]][,s[t+1]] * L[,t])
    
    p.max = which.max(p)
    p[p.max] = max(shift.scale, p[p.max])
    p[-p.max] = p[-p.max]/sum(p[-p.max]+.Machine$double.eps) * (1-p[p.max])

    if(any(is.na(p))) {
      browser()
    }
    
  
    s[t] = sample(x = 1:N, size = 1, prob = p)
    # s[t] = which.max(p)
  }
  
  s
}