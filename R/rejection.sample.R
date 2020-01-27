#' Rejection sampler
#' 
#' @param n number of samples to draw
#' @param log.f log of target density
#' @param log.g log of envelope kernel
#' @param rg function to draw samples from g; should only accept "n" as an 
#'   argument
#' @param a scaling factor for envelope.  a is such that g/a > f everywhere.
#' @param max.it maximum number of tries before failing to generate sample
#' 
#' @return a random draw from f, or NA if sampling fails
#' 
#' @example examples/rejection.sample.R
#' 
rejection.sample = function(n, log.f, log.g, rg, a, max.it = 1e3) {
  
  if(a > 1 | a < 0) {
    stop('Invalid scaling factor; a must lie in (0,1].')
  }
  
  log.a = log(a)
  
  sapply(1:n, function(i) {
    x = NA
    for(j in 1:max.it) {
      y = rg(1)
      if(log(runif(1)) <= log.f(y) - log.g(y) + log.a) {
        x = y
        break
      }
    }
    x
  })
}