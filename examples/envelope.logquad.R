library(truncdist)

# define range for the envelope
xmin = 1
xmax = 10

# define envelope segments
breaks = seq(from = xmin, to = xmax, length.out = 10)

# get number of segments
n = length(breaks)


#
# example 1: f(x) = exp(-x); log f(x) = -x
#

# build envelope
q = envelope.logquad(breaks = breaks, 
                     logf = -breaks[1:(n-1)], 
                     d.logf = rep(-1, n-1), 
                     dd.logf.sup = rep(0, n-1))

# sample from density
samples = q$rquad(n = 1e4)

# compare density functions, and kernel density estimate
curve(dtrunc(x = x, spec = 'exp', a = xmin, b = xmax), from = xmin, to = xmax, 
      ylab = expression(f(x)))
curve(q$dquad(x), from = xmin, to = xmax, col = 2, lty = 2, add = TRUE)
lines(density(samples), col = 3, lty = 2)


#
# example 2: f(x) = exp(sin(x))/C; log f(x) = sin(x) - log(C)
#

# integration constant
C = integrate(function(x) exp(sin(x)), lower = xmin, upper = xmax)$value

# build envelope
q = envelope.logquad(breaks = breaks, 
                     logf = sin(breaks[1:(n-1)]) - log(C), 
                     d.logf = cos(breaks[1:(n-1)]), 
                     dd.logf.sup = sapply(1:(n-1), function(i) {
                       extrema = floor(2*breaks[i + 0:1]/pi) * pi / 2
                       extrema = extrema[
                         (extrema >= breaks[i]) & extrema <= breaks[i+1]
                         ]
                       max(-sin(c(extrema, breaks[i + 0:1])) )
                     }))

# build envelope
anchors = breaks[1:(length(breaks)-1)] + diff(breaks)/2
q2 = envelope.logquad(breaks = breaks, 
                      logf = sin(anchors) - log(C), 
                      d.logf = cos(anchors), 
                      dd.logf.sup = sapply(1:(n-1), function(i) {
                        extrema = floor(2*breaks[i + 0:1]/pi) * pi / 2
                        extrema = extrema[
                          (extrema >= breaks[i]) & extrema <= breaks[i+1]
                          ]
                        max(-sin(c(extrema, breaks[i + 0:1])) )
                      }), 
                      anchors = anchors)

# look at envelope
curve(exp(sin(x))/C, from = xmin, to = xmax, 
      ylab = expression(f(x)))
curve(q$e(x), from = xmin, to = xmax, col = 2, lty = 2, add = TRUE)
curve(q2$e(x), from = xmin, to = xmax, col = 4, lty = 2, add = TRUE)

# sample from density
samples = q2$rquad(n = 1e4)

# compare density to kernel density estimate
curve(q2$dquad(x), from = xmin, to = xmax, ylab = expression(f(x)))
lines(density(samples), col = 3, lty = 2)

# compare CDFs 
curve(sapply(x, function(x) {
  integrate(function(x) exp(sin(x))/C, lower = xmin, upper = x)$value
}), from = xmin, to = xmax, ylab = expression(F(x)))
curve(q2$pquad(x), from = xmin, to = xmax, add = TRUE, col = 2)

# compare density envelopes
curve(exp(sin(x))/C, from = xmin, to = xmax, ylab = expression(f(x)))
curve(q2$dquad(x)*q2$C, from = xmin, to = xmax, add = TRUE, col = 2)

# rejection sampler efficiency when using q2 as an envelope
1/q2$C
