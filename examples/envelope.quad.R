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
