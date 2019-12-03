#' Extract a segment of a fully observed dive trajectory between timepoints
#' 
#' @param depths Complete record of depth bin indices visited
#' @param times Times at which each of \code{depths} was visited
#' @param stages Stages at which each of the \code{depths} was visited
#' @param durations the amount of time spent in each depth bin
#' @param t0 the time to mark the beginning of the segment to extract
#' @param tf the time to mark the end of the segment to extract
#' 
#' #@example examples/dsdive.extract.R
#' 
dsdive.extract = function(depths, times, stages, durations, t0, tf) {
  
  tind = findInterval(c(t0,tf), times, rightmost.closed = TRUE)
  tind.seq = seq(from = tind[1], to = tind[2])
  
  res = list(
    depths = depths[tind.seq],
    times = times[tind.seq],
    stages = stages[tind.seq],
    durations = durations[tind.seq]
  )
  
  class(res) = 'dsobs'
  
  res
}