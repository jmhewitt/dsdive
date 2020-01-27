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

devtools::document()

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
# example 2: f(x) = exp(sin(x)); log f(x) = sin(x)
#


# build envelope
q = envelope.logquad(breaks = breaks, 
                     logf = sin(breaks[1:(n-1)]), 
                     d.logf = cos(breaks[1:(n-1)]), 
                     dd.logf.sup = sapply(1:(n-1), function(i) {
                       max(-sin(c(breaks[i:(i+1)])))
                     }))

# look at envelope
curve(exp(sin(x)), from = xmin, to = xmax, 
      ylab = expression(f(x)))
curve(q$e(x), from = xmin, to = xmax, col = 2, lty = 2, add = TRUE)


# sample from density
samples = q$rquad(n = 1e4)

# compare density to kernel density estimate
curve(q$dquad(x), from = xmin, to = xmax, ylab = expression(f(x)))
lines(density(samples), col = 3, lty = 2)
