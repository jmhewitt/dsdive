#' Initialize a computing environment for use with gibbs sampling
#' 
#'
#' @export
#' 
#' @example examples/makeCompleteSingle.R
#' 
makeCompleteSingle = function(depth.bins, durations, depths, times, stages, 
                              init, priors, t0.dive, model) {
  
  # ensure depth.bins is in a good format
  depth.bins = as.matrix(depth.bins)
  
  # verify input represents a completely or partially observed dive
  if(is.null(durations)) {
    stop(paste('makeCompleteSingle requires complete information about a dive;',
               'Please ensure duration and stage information are provided.'))
  }
  partially.observed = FALSE
  
  # compute initial jacobian
  logJ = params.toList(par = params.toVec(par = init, spec = priors), 
                       spec = priors)$logJ
  
  # set initial log density and trajectory object
  trajectory = list(
    depths = depths,
    stages = stages,
    times = times,
    durations = durations
  )
  
  # extract time at which stage 2 was entered, if any
  stage2.inds = which(trajectory$stages==2)
  if(length(stage2.inds)==0) {
    t.stage2 = NA
  } else {
    t.stage2 = trajectory$times[min(stage2.inds)]
  }
  
  # evaluate initial log-density
  ld = dsdive.ld(depths = trajectory$depths, 
                 durations = trajectory$durations, 
                 times = trajectory$times, stages = trajectory$stages, 
                 beta = init$beta, lambda = init$lambda, 
                 sub.tx = init$sub.tx, surf.tx = init$surf.tx, 
                 depth.bins = depth.bins, t0.dive = t0.dive, 
                 t.stage2 = t.stage2, model = model) + 
    logJ

  # package configuration
  res = list(
    trajectory = trajectory,
    depth.bins = depth.bins, 
    partially.observed = partially.observed,
    ld = ld,
    t.stage2 = t.stage2,
    t0.dive = t0.dive,
    model = model
  )
  
  class(res) = 'dsCompleteSingle'
  
  res
}