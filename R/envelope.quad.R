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
envelope.quad = function(breaks, f, df, ddf.sup, 
                         anchors = breaks[1:(length(breaks)-1)]) {
  
  if(any(f<0)) {
    stop('Probabilistic envelopes can only be built for non-negative fns.')
  }
  
  # build envelope and extract key components
  bound.quad = bound.quad(breaks = breaks, f = f, df = df, ddf.sup = ddf.sup, 
                          anchors = anchors)
  coefs.mat = bound.quad$coefs
  n = nrow(coefs.mat)
  e = bound.quad$e
  
  # build partial integral function, which is a CDF component
  mass.segment = function(i, x) {
    # Integrate the i^th quadratic envelope from start of segment to t = x
    polyval(c(coefs.mat[i,]/3:1, 0), x - anchors[i]) - 
    polyval(c(coefs.mat[i,]/3:1, 0), breaks[i] - anchors[i])
  }
  
  # compute cumulative mass at start of each segment
  mass.cum = c(0, cumsum(sapply(1:n, function(i) mass.segment(i, breaks[i+1]))))
  
  # standardized density function
  dquad = function(x, log = FALSE) {
    r = e(x) / mass.cum[n+1]
    if(log) { log(r) } else { r }
  }
  
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
  
  # inverse CDF associated with envelope
  qquad = function(p) {
    sapply(p, function(p) {
      uniroot(f = function(x) {pquad(x) - p}, lower = breaks[1], 
              upper = breaks[n+1])$root
    })
  }
  
  # inverse-CDF sampling
  rquad = function(n) {
    qquad(runif(n))
  }

 list(pquad = pquad, dquad = dquad, qquad = qquad, rquad = rquad, 
      C = mass.cum[length(mass.cum)])
}