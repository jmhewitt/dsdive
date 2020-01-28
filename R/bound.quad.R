#' Build a piecewise-quadratic envelope for a given function
#' 
#' @param breaks Breakpoints that define the start and end points of each 
#'   interval over which a piecewise-quadratic function will be defined.
#' @param f Value of f at the breakpoints.
#' @param df Derivative of f evaluated at the start of each interval.
#' @param ddf.sup Supremum of second derivative of f over each interval
#' 
#' @example examples/bound.quad.R
#' 
bound.quad = function(breaks, f, df, ddf.sup, 
                      anchors = breaks[1:(length(breaks)-1)]) {
  
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
             polyval(coefs.mat[ind,], x - anchors[ind])
             )
    })
  }

 list(e = e, coefs = coefs.mat, breaks = breaks)
}