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
  
  # specify parameters for the piecewise-quadratic functions
  a = ddf.sup / 2
  b = df
  c = f
  
  coefs.mat = cbind(a, b, c)
  
  # assemble into an envelope function that can be evaluated
  e = function(x) {
    rows = nrow(coefs.mat)
    sapply(x, function(x){
      ind = findInterval(x, breaks)
      ifelse(ind > rows | ind < 1, NA, 
             polyval(coefs.mat[ind,], x - breaks[ind]))
    })
  }
  
 e 
}