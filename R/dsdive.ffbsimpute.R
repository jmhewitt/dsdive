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
#' @example examples/dsdive.ffbsimpute.R
#' @export
#'
dsdive.ffbsimpute = function(depth.bins, depths, times, s0, beta, 
                                        lambda, sub.tx, surf.tx, 
                                        inflation.factor.lambda = 1.1, 
                                        verbose = FALSE, 
                                        t0.dive, t.stages, model, 
                             times.cond = NULL, shift.scale = .9) {
  
  # pre-compute the common max transition rate
  lambda.max = max(outer(lambda, 2 * depth.bins[,2], '/'))
  lambda.thick = lambda.max * inflation.factor.lambda
  
  
  #
  # compute single-stage transition matrices
  #
  
  t1 = dsdive.tx.matrix.simplified(t0 = t0.dive, depth.bins = depth.bins, 
                                   beta = beta, lambda = lambda, 
                                   sub.tx = sub.tx, surf.tx = surf.tx, 
                                   s0 = 1, sf = 1, 
                                   inflation.factor.lambda = 1, 
                                   t0.dive = t0.dive, lambda.max = lambda.thick, 
                                   t.stage2 = NA, model = model)$m
  
  t2 = dsdive.tx.matrix.simplified(t0 = t.stages[1], depth.bins = depth.bins, 
                                   beta = beta, lambda = lambda, 
                                   sub.tx = sub.tx, surf.tx = surf.tx, 
                                   s0 = 2, sf = 2, 
                                   inflation.factor.lambda = 1, 
                                   t0.dive = t0.dive, lambda.max = lambda.thick, 
                                   t.stage2 = t.stages[1], model = model)$m
  
  t3 = dsdive.tx.matrix.simplified(t0 = t.stages[2], depth.bins = depth.bins, 
                                   beta = beta, lambda = lambda, 
                                   sub.tx = sub.tx, surf.tx = surf.tx, 
                                   s0 = 3, sf = 3, 
                                   inflation.factor.lambda = 1, 
                                   t0.dive = t0.dive, lambda.max = lambda.thick, 
                                   t.stage2 = t.stages[1], model = model)$m
  
  
  #
  # sample virtual arrival times from uniformization process
  #
  
  t0 = times[1]
  t.win = times[length(times)] - times[1]
  N = rpois(n = 1, lambda = lambda.thick * t.win)
  W = sort(unique(c(times.cond, t0, t0 + t.win * sort(runif(n = N)))))
  
  
  #
  # build observation vector and transition matrices
  #
    
  N.W = length(W)
  
  W.augmented = c(W,Inf)
  
  O = numeric(N.W)
  Bt = vector('list', N.W)
  
  stages = findInterval(W, t.stages) + 1
  
  for(i in 1:N.W) {
    
    # determine transition matrix
    if(stages[i] == 1) {
      Bt[[i]] = t1
    } else if(stages[i] == 2) {
      Bt[[i]] = t2
    } else {
      Bt[[i]] = t3
    }
    
    # associate observations to time windows
    obs.ind = which((W[i] <= times) & (times < W.augmented[i+1]))
    O[i] = ifelse(length(obs.ind) > 0, depths[obs.ind], NA)
    
    if(!is.na(O[i])) {
      if(O[i] == 1) {
        if(obs.ind > 1) {
          O[i] = nrow(depth.bins) + 1
        }
      }
    }
    
  }
  
  
  #
  # sample and thin transitions
  #
  
  
  # sample transitions
  x = ffbs(Bt = Bt, O = O, shift.scale = shift.scale, W = W, 
           lambda = lambda.thick)
  
  # decode any absorbing returns to the surface
  x[x == (nrow(depth.bins) + 1)] = 1
  
  # remove self-transitions, and include initial state
  true.tx.inds = c(TRUE, diff(x) != 0)
  
  res = list(
    depths = x[true.tx.inds],
    times = W[true.tx.inds],
    stages = stages[true.tx.inds]
  )
  
  res$durations = diff(res$times)
  
  class(res) = 'dsdive'
  
  res
}