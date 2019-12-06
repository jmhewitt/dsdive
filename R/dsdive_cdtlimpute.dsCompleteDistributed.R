#' Conditionally impute trajectories
#'
#' 
dsdive_cdtlimpute.dsCompleteDistributed = function(cfg, params, i) {
  
  # get current log-density
  cfg$ld = dsdive_ld(cfg = cfg, params = params) + params$logJ
  
  # return (updated) cfg
  cfg
}