#' Initialize a computing environment for use with gibbs sampling
#' 
#' 
#'
#' @export
#' 
makeSingle = function(depth.bins, durations, it, depths, times, init, priors,
                      inflation.factor.lambda, t0.dive, verbose = FALSE) {
  
  # ensure depth.bins is in a good format
  depth.bins = as.matrix(depth.bins)
  
  # determine whether input represents a completely or partially observed dive
  partially.observed = is.null(durations)
  
  # compute initial jacobian
  logJ = params.toList(par = params.toVec(par = init, spec = priors), 
                       spec = priors)$logJ
  
  # if necessary, build storage for imputed trajectories, and initialize 
  #   latent trajectory
  if(partially.observed) {
    
    if(verbose) {
      message('Imputing initial trajectory')
    }
    
    # initialize record of latent trajectories
    trace.imputed = vector('list', it)
    
    # sample first latent trajectory
    trajectory = dsdive.fastimpute(M = 1, depth.bins = depth.bins, 
                                   depths = depths, times = times, 
                                   s0 = 1, beta = init$beta, 
                                   lambda = init$lambda, 
                                   sub.tx = init$sub.tx, 
                                   surf.tx = init$surf.tx, 
                                   inflation.factor.lambda = 
                                     inflation.factor.lambda, 
                                   verbose = FALSE, 
                                   precompute.bridges = TRUE, 
                                   t0.dive = t0.dive, 
                                   resample = FALSE)[[1]]
    
    # add log jacobians to sample
    trajectory$ld.true = trajectory$ld.true + logJ
    trajectory$ld = trajectory$ld + logJ
    
    # extract log-density
    ld = trajectory$ld.true
    
    # save sample
    trace.imputed[[1]] = trajectory
    
  } else {
    # otherwise, set initial log density and trajectory object
    trajectory = list(
      depths = depths,
      stages = stages,
      times = times,
      durations = durations
    )
    
    # extract time at which stage 2 was entered, if any
    stage2.inds = which(trajectory$stages==2)
    if(length(stage2.inds)==0) {
      t.stage2tmp = NA
    } else {
      t.stage2tmp = trajectory$times[min(stage2.inds)]
    }
    
    ld = dsdive.ld(depths = trajectory$depths, 
                   durations = trajectory$durations, 
                   times = trajectory$times, stages = trajectory$stages, 
                   beta = init$beta, lambda = init$lambda, 
                   sub.tx = init$sub.tx, surf.tx = init$surf.tx, 
                   depth.bins = depth.bins, t0.dive = t0.dive, 
                   t.stage2 = t.stage2tmp) + 
      logJ
  }
  
  # extract time at which stage 2 was entered, if any
  stage2.inds = which(trajectory$stages==2)
  if(length(stage2.inds)==0) {
    t.stage2 = NA
  } else {
    t.stage2 = trajectory$times[min(stage2.inds)]
  }
  
  #
  # package configuration
  #
  
  res = list(
    trajectory = trajectory,
    depth.bins = depth.bins, 
    partially.observed = partially.observed,
    ld = ld,
    t.stage2 = t.stage2,
    t0.dive = t0.dive
  )
  
  if(partially.observed) {
    res = c(res, list(
      depths = depths,
      times = times,
      trace.imputed = trace.imputed,
      inflation.factor.lambda = inflation.factor.lambda
    ))
  }
  
  class(res) = 'dssingle'
  
  res
}