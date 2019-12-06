#' Evaluate the log-density of a single dive trajectory
#' 
#' Evaluates the log-density of a single dive trajectory when the computing 
#' environment is configured for use on a single core.
#' 
#' @param cfg \code{dssingle} object containing dive trajectory data
#' @param params parameter values at which to evaluate the model's log-density
#'
#' @importFrom snow clusterApply
#' 
dsdive_ld.dsCompleteDistributed = function(cfg, params) {
  
  # distribute likelihood computation
  lds = clusterApply(cl = cfg$cl, x = cfg$ids, fun = function(id) {
    
    # extract local config
    cfg.local = get(x = id, envir = globalenv())
    
    # compute and return log-density
    dsdive_ld(cfg.local)
  })
  
  # return total log-density
  sum(unlist(lds))
}