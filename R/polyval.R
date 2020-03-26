#' Use Horner's method to efficiently evaluate an n degree polynomial
#' 
#' f(x) = a_n * x^n + a_{n-1} * x^{n-1} + ... + a_0
#' 
#' @param coefs polynomial coefficients in decreasing order
#' @param x vector of values of x at which to evaluate polynomial
#' 
#' # @example examples/polyval.R
#' 
polyval = function(coefs, x) {
  
  # extract number of coefficients
  nc = length(coefs) 
  
  # initialize output, y = f(x)
  y = rep(coefs[1], length(x))
  
  # evaluate polynomial in nested sequence
  for(i in 2:nc) {
    y = y*x + coefs[i]
  }
  
  y
}