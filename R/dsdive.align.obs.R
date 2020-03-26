#' Realign observations of a dive to the time scale of the dive
#' 
#' Satellite tag observations of a dive record depths at regular intervals, 
#' irrespective of whether or not an animal is diving.  Given initial and final 
#' offsets \code{offset} and \code{offset.tf}, respectively, this function will 
#' shift the satellite-tag observation \code{times} so that they are referenced 
#' with respect to the start of the dive (i.e., t=0).  Observations in 
#' \code{depths} will be trimmed from the output if they exceed the times 
#' implied by the offsets.
#' 
#' @param depths Indices of observed depth bins
#' @param times Times at which each of \code{depths} was observed
#' @param t.stages Stage transition times for the dive; will be used to compute
#'   the dive stage for each observation
#' @param offset Amount of time by which \code{times} over or under reports the 
#'   true time within the dive of each \code{depths} observation
#' @param offset.tf Amount of time by which the last \code{time} over or under 
#'   reports the true duration of the dive.
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