#' Conditionally impute trajectories
#'
#' 
dsdive_cdtlimpute.dsImputedDistributed = function(cfg, params, i) {
  
  # distribute likelihood computation
  lds = clusterApply(cl = cfg$cl, x = cfg$ids, fun = function(id) {
    
    # extract local config
    cfg.local = get(x = id, envir = globalenv())
    
    # update imputed trajectory
    cfg.local = dsdive_cdtlimpute(cfg = cfg.local, params = params, i = i)
      
    # save changes locally
    assign(x = id, value = cfg.local, envir = globalenv())
    
    # return log-density
    cfg.local$ld
  })
  
  # update current log-density
  cfg$ld = sum(unlist(lds))
  
  # return (updated) cfg
  cfg
}