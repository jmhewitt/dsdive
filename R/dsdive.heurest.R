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
#' @param t0.dive Time at which dive started
#' @param method Method to use in \code{optim} search for initial parameters
#' 
#' # @example examples/heurest.R
#' 
#' @export
#' 
#'
dsdive.heurest = function(depths, times, stages.est, depth.bins, 
                          priors.sd = rep(5, 12), priors.mean = rep(0, 12), 
                          t0.dive, method = 'Nelder-Mead') {
  
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
  times.complete = c(0, cumsum(durations.complete)) + t0.dive
  durations.complete = c(durations.complete, NA)
  stages.complete = c(stages.complete, 3)
  
  # extract time at which stage 2 was entered, if any
  stage2.inds = which(stages.complete==2)
  if(length(stage2.inds)==0) {
    t.stage2tmp = NA
  } else {
    t.stage2tmp = times.complete[min(stage2.inds)]
  }
  
  #
  # initial parameters for model parameters
  #
  
  # estimate diving rate magnitude between observations (m/sec)
  dy = abs(diff(2*depth.bins[depths,1]) / diff(times))
  
  # "moment" estimates for lambda
  stages.tx = which(diff(stages.est) == 1)
  lambda.est = c(mean(dy[1:stages.tx[1]]), 
                 mean(dy[(stages.tx[1]+1):stages.tx[2]]),
                 mean(dy[-(1:stages.tx[2])]))
  
  # ad-hoc estimates for stage transition parameters
  stages.tx.complete = which(diff(stages.complete) == 1)
  t0.dive = times[1]
  tx.mu = log(c(times.complete[stages.tx.complete[1]] - t0.dive, 
                diff(times.complete[stages.tx.complete]))/60)
  tx.win = log((rep(mean(diff(times)) , 2))/60)
  
  #
  # optimize initial parameter estimates
  #
  
  o = optim(par = c(rep(0,3), log(lambda.est), c(tx.mu[1], tx.win[1], 
                                                 tx.mu[2], tx.win[2])), 
            fn = function(theta) {
              theta.list = params.toList(par = theta)
              dsdive.ld(depths = depths.complete, stages = stages.complete, 
                        durations = durations.complete, times = times.complete,
                        beta = theta.list$beta, 
                        lambda = theta.list$lambda, 
                        sub.tx = theta.list$sub.tx,
                        surf.tx = theta.list$surf.tx, depth.bins = depth.bins, 
                        t0.dive = t0.dive, t.stage2 = t.stage2tmp) + 
              sum(dnorm(theta, mean = priors.mean, sd = priors.sd, 
                        log = TRUE)) + theta.list$logJ
            }, method = method, control = list(fnscale = -1))
  
  # package results
  params.toList(par = o$par)
}