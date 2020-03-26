#' Refine a partition of an interval
#' 
#' The input \code{breaks} defines a partition.  \code{refine.partition} will 
#' add evenly spaced points between \code{breaks} such that the maximum 
#' difference between any two points in the output is \code{max.width}.
#' 
#' @param breaks Initial partition of an interval on the real line
#' @param max.width Maximum width of intervals in output partition
#' 
#' # @example examples/refine.partition.R
#' 
refine.partition = function(breaks, max.width) {
  
  n = length(breaks)
  
  do.call(c, lapply(1:n, function(i) {
    if(i==n) {
      breaks[n]
    } else {
      # width of i^th interval
      w = diff(breaks[i + 0:1])
      if(w <= max.width) {
        breaks[i]
      } else {
        # ceiling(w/max.width) == number of internal breakpoints, we add 2 to 
        # n to account for the start and end points
        n = ceiling(w/max.width) + 2
        # but we do not return the endpoint because this is included in the next
        # segment
        seq(from = breaks[i], to = breaks[i+1], length.out = n)[1:(n-1)]
      }
    }
  }))
}