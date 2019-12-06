#' Evaluate the log-density of a single dive trajectory
#' 
#' Evaluates the log-density of a single dive trajectory when the computing 
#' environment is configured for use on a single core.
#' 
#' @param cfg \code{dssingle} object containing dive trajectory data
#' @param params parameter values at which to evaluate the model's log-density
#'
#' 
dsdive_ld.dsImputedSingle = function(cfg, params) {
  
  # extract trajectory
  trajectory = cfg$trajectory
  
  # return log-density
  dsdive.ld(depths = trajectory$depths, durations = trajectory$durations,
            times = trajectory$times, stages = trajectory$stages,
            beta = params$beta, lambda = params$lambda, sub.tx = params$sub.tx,
            surf.tx = params$surf.tx, depth.bins = cfg$depth.bins, 
            t0.dive = cfg$t0.dive, t.stage2 = cfg$t.stage2)
}