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

# compare envelope to target function
curve(exp(x), from = xmin, to = xmax)
curve(q(x), from = xmin, to = xmax, add = TRUE, col = 4, n = 1e3)


#
# example 2: f(x) = exp(-x)
#

# build envelope
q2 = envelope.quad(breaks = breaks, 
                   f = exp(-breaks[1:(n-1)]), 
                   df = -exp(-breaks[1:(n-1)]), 
                   ddf.sup = exp(-breaks[1:(n-1)]))

# compare envelope to target function
curve(exp(-x), from = xmin, to = xmax)
curve(q2(x), from = xmin, to = xmax, add = TRUE, col = 4, n = 1e3)
