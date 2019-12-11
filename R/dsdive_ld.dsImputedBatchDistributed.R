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
dsdive_ld.dsImputedBatchDistributed = function(cfg, params) {
  
  # distribute likelihood computation
  lds = clusterApply(cl = cfg$cl, x = cfg$cluster.ids, fun = function(id) {
    
    # extract local config
    cfg.local = get(x = id, envir = globalenv())
    
    # compute and return log-density
    if(!is.null(cfg.local)) {
      dsdive_ld(cfg = cfg.local, params = params)
    } else {
      0
    }
    
  })
  
  # return total log-density
  sum(unlist(lds))
}