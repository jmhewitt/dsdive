#' Heuristic estimate of model parameters given observations
#' 
#' Builds heuristic estimates of dive model parameters, which is useful for 
#' initializing methods that optimize or explore parameter surfaces and 
#' posterior distributions.  
#' 
#' The heuristic estimate begins by linearly imputing a dive trajectory that 
#' is consistent with observations, then optimizes model parameters from the 
#' basic imputed trajectory.  Note that the imputed trajectory is likely to be 
#' very unrealistic, especially for datasets with sparse temporal observations.
#'
#' @param depths Record of observed depth bin indices
#' @param times Times at which depth bin observations were made
#' @param stages.est Vector of guesses for which dive stage the trajectory was 
#'   in at each observation
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param priors Standard deviations for mean-0 normal priors on model 
#'   parameters, after transformation to real line
#' @param sub.tx1 Index of minimum depth bin at which transition to intermediate 
#'   behaviors stage can occur.
#' @param t0.dive Time at which dive started
#' 
#' # @example examples/heurest.R
#' 
#' @export
#' 
#'
dsdive.heurest = function(depths, times, stages.est, depth.bins, 
                          priors = rep(5, 12), sub.tx1, t0.dive) {
  
  #
  # linearly impute a complete trajectory
  #
  
  nt = length(times)
  
  # impute direct paths, equal time bins, and stages between observations
  depths.complete = depths[1]
  durations.complete = c()
  stages.complete = c()
  for(i in 2:nt) {
    # impute direct path between observations
    path = (depths.complete[length(depths.complete)]:depths[i])[-1]
    # associate equal time bins for the path
    np = length(path)
    path.times = rep(diff(times[(i-1):i]) / np, np)
    # assign stages to path
    path.stages = rep(stages.est[i-1], np)
    # record imputation
    depths.complete = c(depths.complete, path)
    durations.complete = c(durations.complete, path.times)
    stages.complete = c(stages.complete, path.stages)
  }
  
  # finish assembling complete trajectory
  times.complete = c(0, cumsum(durations.complete))
  durations.complete = c(durations.complete, NA)
  stages.complete = c(stages.complete, 3)
  
  
  #
  # initial parameters for model parameters
  #
  
  # estimate diving rate magnitude between observations (m/sec)
  dy = abs(diff(2*depth.bins[depths,2]) / diff(times))
  
  # "moment" estimates for lambda
  stages.tx = which(diff(stages.est) == 1)
  lambda.est = c(mean(dy[1:stages.tx[1]]), 
                 mean(dy[(stages.tx[1]+1):stages.tx[2]]),
                 mean(dy[-(1:stages.tx[2])]))
  
  # "moment" estimate for the probability of transitioning from primary dive
  sub.tx2 = 1/(
    min(which(stages.complete == 2)) - min(which(depths.complete >= sub.tx1))
  )
  
  
  #
  # optimize initial parameter estimates
  #
  
  o = optim(par = c(rep(0,6), log(lambda.est), qlogis(sub.tx2), rep(0,2)), 
            fn = function(theta) {
              dsdive.ld(depths = depths.complete, stages = stages.complete, 
                        durations = durations.complete, times = times.complete,
                        beta = matrix(theta[1:6], nrow = 2), 
                        lambda = exp(theta[7:9]), 
                        sub.tx = c(sub.tx1, plogis(theta[10])),
                        surf.tx = theta[11:12], depth.bins = depth.bins, 
                        t0.dive = t0.dive) + 
              sum(dnorm(theta, sd = priors, log = TRUE))
            }, method = 'BFGS', control = list(fnscale = -1))
  
  # package results
  list(
    beta = matrix(o$par[1:6], nrow = 2),
    lambda = exp(o$par[7:9]),
    sub.tx = c(sub.tx1, plogis(o$par[10])),
    surf.tx = o$par[11:12]
  )
}