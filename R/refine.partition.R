#' Refine a partition of an interval
#' 
#' The input \code{breaks} defines a partition.  \code{refine.partition} will 
#' add evenly spaced points between \code{breaks} such that the maximum 
#' difference between any two points in the output is \code{max.width}.
#' 
#' @example examples/refine.partition.R
#' 
refine.partition = function(breaks, max.width) {
  
  c(do.call(c, sapply(1:(length(breaks)-1), function(i) {
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
  })), breaks[length(breaks)])
  
}