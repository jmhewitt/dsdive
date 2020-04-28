#' Approximate slope and curvature for each polynomial segment of an envelope
#' 
#' @param breaks Breakpoints that define the start and end points of each 
#'   interval over which a piecewise-quadratic function will be defined.
#' @param lp.breaks Value of log(f) at the breakpoints.
#' @param lp.anchors Value of log(f) at the anchors.
#' @param anchors Locations within the intervals (defined by \code{breaks}) at 
#'   which the polynomial approximation will be centered
#' 
# @example examples/envelope.approx.R
#' 
envelope.approx = function(breaks, anchors, lp.breaks, lp.anchors) {
  
  envelope = sapply(1:length(anchors), function(i) {
    # extract time-coordinates for interval
    L.x = breaks[i]
    M.x = anchors[i]
    U.x = breaks[i+1]
    # extract log-posterior for interval
    L = lp.breaks[i]
    M = lp.anchors[i]
    U = lp.breaks[i+1]
    # the in/finiteness of the points determines processing
    L.undef = !is.finite(L)
    M.undef = !is.finite(M)
    U.undef = !is.finite(U)
    
    # assume lp = -Inf throughout entire interval
    if(L.undef & M.undef & U.undef) {
      c(0, 0, -Inf)
    }
    # assume lp = -Inf at start and midpoint, but ends finite
    else if(L.undef & M.undef) {
      c(0, 0, U)
    }
    # assume lp is only -Inf at start point
    else if(L.undef) {
      c(0, (U-M)/(U.x-M.x), M)
    }
    # assume lp is -Inf at end and midpoint, but starts finite
    else if(U.undef & M.undef) {
      c(0, 0, L)
    }
    # assume lp is only -Inf at end point
    else if(U.undef) {
      c(0, (L-M)/(L.x-M.x), M)
    }
    # assume lp is finite in entire interval
    else {
      b = (M-L)/(M.x-L.x)
      x = U.x - M.x
      # a = (U - M - b * x)/x^2
      a = 0
      c(2*a, b, M)
    }
  })
  
  list(
    logf = envelope[3,],
    d.logf = envelope[2,],
    dd.logf.sup = envelope[1,]
  )
}