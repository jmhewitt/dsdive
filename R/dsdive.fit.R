#' Likelihood for completely observed dive trajectories
#'
#' @param depths Depth bin indices in which the trajectory was observed
#' @param durations Amount of time spent in each depth bin
#' @param times Times at which each depth bin was entered
#' @param stages Record of dive stages at each depth bin
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
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param t0.dive Time at which dive started
#' @param d0.last If the depth bin that proceeded the first depth bin in 
#'   \code{depths}.  If the trajectory to be analyzed was started at the 
#'   surface, then set \code{c0.last=NULL}.
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#' 
#' @example examples/ld.R
#' 
#' @importFrom stats dexp dbinom
#' 
#' @export
#' 
dsdive.fit = function(
  dives.obs, cl, beta.init, lambda.init, verbose = FALSE, T1.prior.params, 
  T2.prior.params, pi1.prior, pi2.prior, lambda1.prior, lambda2.prior, 
  lambda3.prior, max.width, it, inflation.factor.lambda = 1.1, 
  checkpoint.interval = Inf, checkpoint.function = function(x, ...) {}, 
  method.N = 'exact', N.max = 100, crash.function = function(x, ...) {}, 
  truncpois.scale = .8) {
  
  impute = function(depth.bins, depths, times, beta, lambda, t.stages) {
    dsdive.impute_segments(
      depth.bins = depth.bins, depths = depths, times = times, beta = beta, 
      lambda = lambda, s0 = 1, verbose = FALSE, method.N = method.N,
      inflation.factor.lambda = inflation.factor.lambda, N.max = N.max, 
      t.sbreaks = t.stages, truncpois.scale = truncpois.scale)
  }
  
  impute.gibbs = function(depth.bins, depths, times, beta, lambda, 
                          imputed.cond) {
    impute(depth.bins, depths, times, beta, lambda, imputed.cond$t.stages)
  }
  
  res = dsdive.gibbs(dives.obs = dives.obs, cl = cl, impute.init = impute, 
                     impute.gibbs = impute.gibbs,
                     init = list(beta = beta.init, lambda = lambda.init), 
                     verbose = verbose, maxit = it, 
                     T1.prior.params = T1.prior.params, 
                     T2.prior.params = T2.prior.params, pi1.prior = pi1.prior, 
                     pi2.prior = pi2.prior, lambda1.prior = lambda1.prior, 
                     lambda2.prior = lambda2.prior, 
                     lambda3.prior = lambda3.prior, 
                     checkpoint.fn = checkpoint.function, 
                     checkpoint.interval = checkpoint.interval, 
                     crash.fn = crash.function, max.width = max.width)
   
  res 
}