#' Simulate observing a fully observed dive trajectory at specified timepoints
#' 
#' @example examples/txparams.R
#' 
#' @export
#' 
dsdive.observe = function(depths, times, t.obs) {
  list(
    depths = depths[findInterval(t.obs, times)],
    times = t.obs
  )
}