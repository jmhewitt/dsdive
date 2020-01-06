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
#' @export
#'
dsdive.batchimpute = function(depth.bins, depths, times, s0, beta, 
                                        lambda, sub.tx, surf.tx, tx.max,
                                        inflation.factor.lambda = 1.1, 
                                        verbose = FALSE, 
                                        t0.dive, t.stages, model) {
  
  # pre-compute the common max transition rate
  lambda.max = max(outer(lambda, 2 * depth.bins[,2], '/'))
  
  #
  # impute segments
  #
  
  # find stage breaks
  s1.obs = times < t.stages[1]
  s2.obs = (times < t.stages[2]) & (!s1.obs) 
  s3.obs = times > t.stages[2]
  
  # convert to indices
  s1.inds = which(s1.obs)
  s2.inds = which(s2.obs)
  s3.inds = which(s3.obs)
  
  # impute stage 1 - stage 1 transitions
  s1.imputed = dsdive.batchedbridge(depth.bins = depth.bins, 
                                    depths = depths[s1.inds], d0.last = NULL, 
                                    times = times[s1.inds], tx.max = tx.max, 
                                    beta = beta, lambda = lambda, 
                                    sub.tx = sub.tx, surf.tx = surf.tx, 
                                    s0 = 1, sf = 1, 
                                    inflation.factor.lambda = 
                                      inflation.factor.lambda, verbose = FALSE, 
                                    lambda.max = lambda.max, t0.dive = t0.dive, 
                                    t.stages = t.stages, model = model)
  
  # impute stage 1 - stage 2 transitions
  s12.inds = s1.inds[length(s1.inds)] + 0:1
  s12.imputed = dsdive.batchedbridge(depth.bins = depth.bins, 
                                    depths = depths[s12.inds], 
                                    d0.last = s1.imputed$depths[
                                      length(s1.imputed$depths)-1], 
                                    times = times[s12.inds], tx.max = tx.max, 
                                    beta = beta, lambda = lambda, 
                                    sub.tx = sub.tx, surf.tx = surf.tx, 
                                    s0 = 1, sf = 2, 
                                    inflation.factor.lambda = 
                                      inflation.factor.lambda, verbose = FALSE, 
                                    lambda.max = lambda.max, t0.dive = t0.dive, 
                                    t.stages = t.stages, model = model)
  
  # impute stage 2 - stage 2 transitions
  s2.imputed = dsdive.batchedbridge(depth.bins = depth.bins, 
                                    depths = depths[s2.inds], 
                                    d0.last = s12.imputed$depths[
                                      length(s12.imputed$depths)-1], 
                                     times = times[s2.inds], tx.max = tx.max, 
                                     beta = beta, lambda = lambda, 
                                     sub.tx = sub.tx, surf.tx = surf.tx, 
                                     s0 = 2, sf = 2, 
                                     inflation.factor.lambda = 
                                       inflation.factor.lambda, verbose = FALSE, 
                                     lambda.max = lambda.max, t0.dive = t0.dive, 
                                     t.stages = t.stages, model = model)
  
  # impute stage 2 - stage 3 transitions
  s23.inds = s2.inds[length(s2.inds)] + 0:1
  s23.imputed = dsdive.batchedbridge(depth.bins = depth.bins, 
                                    depths = depths[s23.inds], 
                                    d0.last = s2.imputed$depths[
                                      length(s2.imputed$depths)-1], 
                                    times = times[s23.inds], tx.max = tx.max, 
                                    beta = beta, lambda = lambda, 
                                    sub.tx = sub.tx, surf.tx = surf.tx, 
                                    s0 = 2, sf = 3, 
                                    inflation.factor.lambda = 
                                      inflation.factor.lambda, verbose = FALSE, 
                                    lambda.max = lambda.max, t0.dive = t0.dive, 
                                    t.stages = t.stages, model = model)
  
  # impute stage 3 - stage 3 transitions
  s3.imputed = dsdive.batchedbridge(depth.bins = depth.bins, 
                                    depths = depths[s3.inds], 
                                    d0.last = s23.imputed$depths[
                                      length(s23.imputed$depths)-1], 
                                    times = times[s3.inds], tx.max = tx.max, 
                                    beta = beta, lambda = lambda, 
                                    sub.tx = sub.tx, surf.tx = surf.tx, 
                                    s0 = 3, sf = 3, 
                                    inflation.factor.lambda = 
                                      inflation.factor.lambda, verbose = FALSE, 
                                    lambda.max = lambda.max, t0.dive = t0.dive, 
                                    t.stages = t.stages, model = model)
  
  #
  # merge segments
  #
  
  res = list(
    depths = c(s1.imputed$depths, 
               s12.imputed$depths[-1],
               s2.imputed$depths[-1],
               s23.imputed$depths[-1],
               s3.imputed$depths[-1]),
    stages = c(s1.imputed$stages, 
               s12.imputed$stages[-1],
               s2.imputed$stages[-1],
               s23.imputed$stages[-1],
               s3.imputed$stages[-1]),
    times = c(s1.imputed$times, 
               s12.imputed$times[-1],
               s2.imputed$times[-1],
               s23.imputed$times[-1],
               s3.imputed$times[-1])
  )
  
  res$durations = diff(res$times)
  
  class(res) = 'dsdive'
  
  res
}