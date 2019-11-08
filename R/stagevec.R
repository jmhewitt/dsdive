#' Convenience function for building a vector for dive stages
#' 
#' Several functions in the \code{dsdive} package have a record of dive 
#' stages as one of the function's arguments.  However, this information is not 
#' always available.  For example, dive stage information is a latent variable
#' so is not included with observational data.  
#' 
#' The \code{stagevec} function makes it easy to generate a vector of dive 
#' stages, like \code{c(1,1,1,2,2,2,2,3,3)}.  Users only need to provide a 
#' target length \code{length.out} and the indices at which the stage 
#' changes \code{breaks.}
#' 
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