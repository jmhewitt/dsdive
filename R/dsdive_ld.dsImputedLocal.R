#' Evaluate the log-density of a single dive trajectory
#' 
#' Evaluates the log-density of a single dive trajectory when the computing 
#' environment is configured for use on a single core.
#' 
#' @param cfg \code{dssingle} object containing dive trajectory data
#' @param params parameter values at which to evaluate the model's log-density
#'
#' 
dsdive_ld.dsImputedLocal = function(cfg, params) {
  
  # distribute likelihood computation
  lds = lapply(cfg$ids, function(id) {
    
    # extract local config
    cfg.local = get(x = id, envir = globalenv())
    
    # compute and return log-density
    dsdive_ld(cfg = cfg.local, params = params)
  })
  
  # return total log-density
  sum(unlist(lds))
}