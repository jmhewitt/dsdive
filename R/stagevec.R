#' Convenience function for building a vector for dive stages
#' 
#' @param length.out Length of vector to generate
#' @param breaks Location of changepoints
#' 
#' @examples
#' 
#' stagevec(length.out = 30, breaks = c(7, 20))
#' 
#' @export
#' 
stagevec = function(length.out, breaks) {
  rep(1:3, c(breaks[1] - 1, breaks[2] - breaks[1], length.out + 1 - breaks[2]) )
}