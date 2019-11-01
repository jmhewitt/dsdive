#' Simulate observing a fully observed dive trajectory at specified timepoints
#' 
#' @param depths Complete record of depth bins visited
#' @param times Times at which each of \code{depths} was visited
#' @param t.obs Times at which the trajectory should be observed
#' 
#' @example examples/observe.R
#' 
#' @export
#' 
dsdive.observe = function(depths, times, t.obs) {
  list(
    depths = depths[findInterval(t.obs, times)],
    times = t.obs
  )
}