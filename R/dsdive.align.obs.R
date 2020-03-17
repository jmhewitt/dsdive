#' Extract a segment of a fully observed dive trajectory between timepoints
#' 
#' @param depths Complete record of depth bin indices visited
#' @param times Times at which each of \code{depths} was visited
#' @param stages Stages at which each of the \code{depths} was visited
#' @param durations the amount of time spent in each depth bin
#' @param t0 the time to mark the beginning of the segment to extract
#' @param tf the time to mark the end of the segment to extract
#' 
#' @example examples/dsdive.align.obs.R
#' 
#' @export
#' 
dsdive.align.obs = function(depths, times, t.stages, offset, offset.tf) {
  
  # realign observation times with beginning of dive
  times.aligned = times - offset
  
  # remove observations from before the dive begins
  nonneg.times = times.aligned >= 0
  times.aligned = times.aligned[nonneg.times]
  depths.aligned = depths[nonneg.times]
  
  # dives start at the surface; add this if necessary
  if(times.aligned[1] > 0) {
    times.aligned = c(0, times.aligned)
    depths.aligned = c(1, depths.aligned)
  }
  
  # identify observations at or after when the dive ends
  tf = times.aligned[length(times.aligned)] - offset.tf
  times.within = times.aligned < tf
  
  # remove observations after the dive end and add final surface obs
  times.aligned = c(times.aligned[times.within], tf)
  depths.aligned = c(depths.aligned[times.within], 1)
  
  #
  # package dive
  #
  
  aligned = list(
    times = times.aligned,
    depths = depths.aligned,
    stages = findInterval(times, t.stages) + 1
  )
  
  class(aligned) = 'dsobs'
  
  aligned
}