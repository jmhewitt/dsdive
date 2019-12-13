#' Conditionally impute trajectories
#'
#' 
dsdive_cdtlimpute.dsImputedSingle = function(cfg, params, i) {
  
  # extract current times of stage break points
  stage.times = cfg$trajectory$times[which(diff(cfg$trajectory$stages)==1)]
  
  # propose trajectory
  prop = dsdive.fastimpute(M = 2, depth.bins = cfg$depth.bins, 
                           depths = cfg$depths, times = cfg$times, s0 = 1, 
                           beta = params$beta, lambda = params$lambda, 
                           sub.tx = params$sub.tx, surf.tx = params$surf.tx, 
                           inflation.factor.lambda = 
                             cfg$inflation.factor.lambda, verbose = FALSE, 
                           precompute.bridges = TRUE, t0.dive = cfg$t0.dive, 
                           trajectory.conditional = cfg$trajectory,
                           model = cfg$model)
  
  # get raw sampling weights
  W = exp(sapply(prop, function(p) p$w))
  # correct for edge cases in weights
  W[is.infinite(W)] = sign(W[is.infinite(W)])
  W[W==-1] = 0
  
  # sample dive 
  cfg$trajectory = prop[[sample(x = 2, size = 1, prob = W)]]
  
  # sample stage transition times
  if(cfg$stages.conditional) {
    
    stage.breaks =  c(
      max(which(cfg$trajectory$times <= stage.times[1])),
      max(which(cfg$trajectory$times <= stage.times[2]))
    )
      
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
  
  # extract time at which stage 2 was entered, if any
  stage2.inds = which(cfg$trajectory$stages==2)
  if(length(stage2.inds)==0) {
    cfg$t.stage2 = NA
  } else {
    cfg$t.stage2 = cfg$trajectory$times[min(stage2.inds)]
  }
  
  # add jacobians
  cfg$trajectory$ld.true = cfg$trajectory$ld.true + params$logJ
  cfg$trajectory$ld = cfg$trajectory$ld + params$logJ
  
  # update log-density of proposal and trace
  cfg$trace.imputed[i] = list(cfg$trajectory)
  cfg$ld = dsdive.ld(depths = cfg$trajectory$depths, 
                     durations = cfg$trajectory$durations, 
                     times = cfg$trajectory$times, 
                     stages = cfg$trajectory$stages, beta = params$beta, 
                     lambda = params$lambda, sub.tx = params$sub.tx, 
                     surf.tx = params$surf.tx, depth.bins = cfg$depth.bins, 
                     t0.dive = cfg$t0.dive, t.stage2 = cfg$t.stage2, 
                     model = cfg$model) + 
    params$logJ
  
  # return (updated) cfg
  cfg
}