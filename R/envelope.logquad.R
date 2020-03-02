#' Build a piecewise-quadratic envelope given log of a function
#' 
#' @param breaks Breakpoints that define the start and end points of each 
#'   interval over which a piecewise-quadratic function will be defined.
#' @param f Value of f at the breakpoints.
#' @param df Derivative of f evaluated at the start of each interval.
#' @param ddf.sup Supremum of second derivative of f over each interval
#' 
#' @importFrom pracma erfi
#' 
#' @example examples/envelope.logquad.R
#' 
envelope.logquad = function(breaks, logf, d.logf, dd.logf.sup,
                            anchors = breaks[1:(length(breaks)-1)]) {
  
  # build envelope for log(f) and extract key components
  bound.logquad = bound.quad(breaks = breaks, f = logf, df = d.logf, 
                     ddf.sup = dd.logf.sup, anchors = anchors)
  coefs.mat = bound.logquad$coefs
  e.log = bound.logquad$e
  n = nrow(coefs.mat)

  # build partial integral function, which is a CDF component
  mass.segment = function(i, x) {
    # Integrate the i^th quadratic envelope from start of segment to t = x
    
    # extract quadratic coefficients
    a = coefs.mat[i,1] + 0i
    b = coefs.mat[i,2]
    c = coefs.mat[i,3]
    
    # account for mass when log-fn is -Inf
    if(is.infinite(c)) {
      if(c<0) {
        if(all(Re(a) < Inf, b < Inf)) {
          0
        }
      } else {
        NaN
      }
    } 
    # account for mass when log-fn is finite
    else {
      if(a==0) {
        exp(c) / b * ( exp(b * (x - anchors[i])) - 
                         exp(b * (breaks[i] - anchors[i])) )
      } else {
        Re(sqrt(pi) * exp(c-b^2/(4*a)) / 2 / sqrt(a) * (
          erfi((2*a*(x-anchors[i]) + b)/(2*sqrt(a))) -
            erfi((2*a*(breaks[i]-anchors[i]) + b)/(2*sqrt(a)))
        ))
      }
    }
  }
    
    

  # compute cumulative mass at start of each segment
  mass.cum = c(0, cumsum(sapply(1:n, function(i) mass.segment(i, breaks[i+1]))))
  
  # standardized density function
  dquad = function(x, log = FALSE) {
    r = exp(e.log(x)) / mass.cum[n+1]
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
      e = function(x) exp(e.log(x)), e.log = e.log, 
      C = mass.cum[n+1], mass.cum = mass.cum)
}