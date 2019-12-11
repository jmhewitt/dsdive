#' Conditionally impute trajectories
#'
#' 
dsdive_cdtlimpute.dsImputedLocal = function(cfg, params, i) {
  
  # distribute update
  lapply(cfg$ids, function(id) {
    
    # extract local config
    cfg.local = get(x = id, envir = globalenv())
    
    # update imputed trajectory
    cfg.local = dsdive_cdtlimpute(cfg = cfg.local, params = params, i = i)
      
    # save changes locally
    assign(x = id, value = cfg.local, envir = globalenv())
  })
  
  # update current log-density
  cfg$ld = dsdive_ld(cfg = cfg, params = params) + params$logJ
  
  # return (updated) cfg
  cfg
}