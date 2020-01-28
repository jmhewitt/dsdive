# define range for the envelope
xmin = 1
xmax = 10

# define envelope segments
breaks = seq(from = xmin, to = xmax, length.out = 10)

# get number of segments
n = length(breaks)


#
# example 1: f(x) = exp(x)
#

# build envelope
q = envelope.quad(breaks = breaks, 
                  f = exp(breaks[1:(n-1)]), 
                  df = exp(breaks[1:(n-1)]), 
                  ddf.sup = exp(breaks[-1]))

# look at equivalent CDF over the region
curve(q$pquad(x), from = xmin - 1, to = xmax + 1, ylab = expression(F(x)))
abline(v = c(xmin, xmax), lty = 3)

# look at inverse CDF 
curve(q$qquad(x), ylab = expression(F^-1*(x)))

# sample from density
samples = q$rquad(n = 1e4)

# look at standardized density, compare with kernel density estimate
curve(q$dquad(x), from = xmin, to = xmax, ylab = expression(f(x)))
lines(density(samples), col = 2)


#
# example 2: f(x) = (sin(x) + 1) / C
#

# integration constant
C = xmax - xmin + cos(xmin) - cos(xmax)

# build envelope with shifted anchors
anchors = breaks[1:(length(breaks)-1)] + diff(breaks)/2
q2 = envelope.quad(breaks = breaks, 
                   f = (sin(anchors) + 1)/C, 
                   df = cos(anchors)/C, 
                   ddf.sup = sapply(1:(n-1), function(i) {
                    extrema = floor(2*breaks[i + 0:1]/pi) * pi / 2
                    extrema = extrema[
                      (extrema >= breaks[i]) & extrema <= breaks[i+1]
                      ]
                    max(-sin(c(extrema, breaks[i + 0:1]))/C )
                  }), 
                  anchors = anchors)

# compare envelope to target function
curve((sin(x) + 1)/C, from = xmin, to = xmax, n = 1e3)
curve(q2$dquad(x)*q2$C, from = xmin, to = xmax, add = TRUE, col = 2, n = 1e3)

# efficiency of a rejection sampler that uses the envelope q2 
1/q2$C
