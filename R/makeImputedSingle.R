#' Initialize a computing environment for use with gibbs sampling
#' 
#'
#' @param stages.conditional If \code{FALSE}, then dive trajectories and stage
#'   transition times will be jointly sampled, otherwise the stage transition 
#'   times will be sampled conditionally given the complete sequence of depth 
#'   bins and durations.
#' 
#' @export
#' 
#' @example examples/makeImputedSingle.R
#' 
makeImputedSingle = function(depth.bins, it, depths, times, init, 
                             priors, inflation.factor.lambda, t0.dive, 
                             verbose = FALSE, model, 
                             stages.conditional) {
  
  # ensure depth.bins is in a good format
  depth.bins = as.matrix(depth.bins)
  
  # set partial observation
  partially.observed = TRUE
  
  # compute initial jacobian
  logJ = params.toList(par = params.toVec(par = init, spec = priors), 
                       spec = priors)$logJ
  
  if(verbose) {
    message('Imputing initial trajectory')
  }
  
  # initialize record of latent trajectories
  trace.imputed = vector('list', it + 1)
  
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
                                 resample = FALSE, model = model)[[1]]
  
  # sample stage transition times
  if(cfg$stages.conditional) {
    
    stage.breaks =  c(1, round(length(cfg$trajectory$stages)/2))
    
    dens = dsdive.ld.stages(breaks = stage.breaks, fixed.ind = 2, 
                            beta = params$beta, lambda = params$lambda, 
                            sub.tx = params$sub.tx, surf.tx = params$surf.tx, 
                            depths = cfg$trajectory$depths, 
                            durations = cfg$trajectory$durations, 
                            times = cfg$trajectory$times, 
                            depth.bins = cfg$depth.bins, t0.dive = cfg$t0.dive, 
                            t.stage2 = cfg$t.stage2, model = cfg$model)
    stage.breaks[1] = sample(x = dens$x, size = 1, prob = dens$prob)
    
    dens = dsdive.ld.stages(breaks = stage.breaks, fixed.ind = 1, 
                            beta = params$beta, lambda = params$lambda, 
                            sub.tx = params$sub.tx, surf.tx = params$surf.tx, 
                            depths = cfg$trajectory$depths, 
                            durations = cfg$trajectory$durations, 
                            times = cfg$trajectory$times, 
                            depth.bins = cfg$depth.bins, t0.dive = cfg$t0.dive, 
                            t.stage2 = cfg$t.stage2, model = cfg$model)
    stage.breaks[2] = sample(x = dens$x, size = 1, prob = dens$prob)
    
    # update stage vector
    cfg$trajectory$stages = stagevec(length.out = length(cfg$trajectory$stages), 
                                     breaks = stage.breaks)
  }
  
  
  # add log jacobians to sample
  trajectory$ld.true = trajectory$ld.true + logJ
  trajectory$ld = trajectory$ld + logJ
  
  # extract log-density
  ld = trajectory$ld.true
  
  # save sample
  trace.imputed[[1]] = trajectory
  
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
    t0.dive = t0.dive,
    depths = depths,
    times = times,
    trace.imputed = trace.imputed,
    inflation.factor.lambda = inflation.factor.lambda,
    model = model,
    stages.conditional = stages.conditional
  )
  
  class(res) = 'dsImputedSingle'
  
  res
}