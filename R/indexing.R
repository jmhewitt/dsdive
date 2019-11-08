#' Converts (x,y,z) coordinates to an array index
#' 
#' Used for getting the flattened index for an array with dimensions
#' \eqn{x.max x y.max x n}, where \eqn{n} is arbitrary.
#' 
#' @param x x coordinate, ranges between 1 and \code{x.max}
#' @param y y coordinate, ranges between 1 and \code{y.max}
#' @param z z coordinate, must be at least 1
#' @param x.max maximum value of x coordinates
#' @param y.max maximum value of y coordinates
#' 
#' @examples
#' toInd(x = 60, y = 13, z = 1, x.max = 100, y.max = 100)
#' 
#' @export
#' 
toInd = function(x, y, z, x.max, y.max) {
  x + x.max * (y-1 + y.max * (z-1))
}

#' Converts array indices to (x,y,z) coordinates
#' 
#' Used for converting a flattened index to the \eqn{(x,y,z)} indices for an 
#' array with dimensions \eqn{x.max x y.max x n}, where \eqn{n} is arbitrary.
#' 
#' @param ind index in an array
#' @param x.max maximum value of x coordinates
#' @param y.max maximum value of y coordinates
#' 
#' @examples
#' fromInd(ind = 100, x.max = 100, y.max = 100)
#' 
#' @export
#' 
fromInd = function(ind, x.max, y.max) {
  z = ind / x.max / y.max
  z = ifelse(z == round(z), z, floor(z) + 1)
  y = (ind - x.max * y.max * (z-1)) / x.max
  y = ifelse(y == round(y), y, floor(y) + 1)
  x = ind - x.max * (y-1 + y.max * (z-1))
  c(x,y,z)
}