#' Build a piecewise-quadratic envelope for a given function
#' 
#' @param breaks Breakpoints that define the start and end points of each 
#'   interval over which a piecewise-quadratic function will be defined.
#' @param f Value of f at the breakpoints.
#' @param df Derivative of f evaluated at the start of each interval.
#' @param ddf.sup Supremum of second derivative of f over each interval
#' 
#' @example examples/envelope.quad.R
#' 
envelope.quad = function(breaks, f, df, ddf.sup) {
  
  # number of segments for piecewise-quadratic envelope
  n = length(f)
  
  # specify parameters for the piecewise-quadratic functions
  a = ddf.sup / 2
  b = df
  c = f
  
  # assemble parameters into matrix, for simpler access
  coefs.mat = cbind(a, b, c)
  
  # assemble an envelope function, which is an unnormalized density
  e = function(x) {
    sapply(x, function(x){
      ind = findInterval(x, breaks)
      ifelse(ind > n | ind < 1, NA, 
             polyval(coefs.mat[ind,], x - breaks[ind])
             )
    })
  }
  
  # build partial integral function, which is a CDF component
  mass.segment = function(i, x) {
    # Integrate the i^th quadratic envelope from start of segment to t = x
    polyval(c(coefs.mat[i,]/3:1, 0), x - breaks[i])
  }
  
  # compute cumulative mass at start of each segment
  mass.cum = c(0, cumsum(sapply(1:n, function(i) mass.segment(i, breaks[i+1]))))
  
  # CDF associated with piecewise-quadratic envelope (defined for all x \in R)
  pquad = function(x, log = FALSE) {
    r = sapply(x, function(x) {
      ind = findInterval(x, breaks, rightmost.closed = TRUE)
      if(ind == 0) { 0 }
      else if(ind == n + 1) { 1 } 
      else { (mass.cum[ind] + mass.segment(ind, x)) / mass.cum[n+1] }
    })
    if(log) { log(r) } else { r }
  }

 list(e = e, pquad = pquad)
}