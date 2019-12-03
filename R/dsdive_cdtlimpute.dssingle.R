#' Conditionally impute trajectories
#'
#' 
dsdive_cdtlimpute.dssingle = function(cfg, params, i) {
  
  # only update dive trajectories if they are partially observed
  if(cfg$partially.observed) {
    # propose trajectory
    prop = dsdive.fastimpute(M = 2, depth.bins = cfg$depth.bins, 
                             depths = cfg$depths, times = cfg$times, s0 = 1, 
                             beta = params$beta, lambda = params$lambda, 
                             sub.tx = params$sub.tx, surf.tx = params$surf.tx, 
                             inflation.factor.lambda = 
                               cfg$inflation.factor.lambda, verbose = FALSE, 
                             precompute.bridges = TRUE, t0.dive = cfg$t0.dive, 
                             trajectory.conditional = cfg$trajectory)
    
    # get raw sampling weights
    W = exp(sapply(prop, function(p) p$w))
    # correct for edge cases in weights
    W[is.infinite(W)] = sign(W[is.infinite(W)])
    W[W==-1] = 0
    
    # sample dive 
    cfg$trajectory = prop[[sample(x = 2, size = 1, prob = W)]]
    
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
    cfg$ld = cfg$trajectory$ld.true
  }
  
  # return (updated) cfg
  cfg
}