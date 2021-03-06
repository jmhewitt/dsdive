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
q = bound.quad(breaks = breaks, 
               f = exp(breaks[1:(n-1)]), 
               df = exp(breaks[1:(n-1)]), 
               ddf.sup = exp(breaks[-1]))

# compare envelope to target function
curve(exp(x), from = xmin, to = xmax)
curve(q$e(x), from = xmin, to = xmax, add = TRUE, col = 4, n = 1e3)


#
# example 2: f(x) = exp(-x)
#

# build envelope
q2 = bound.quad(breaks = breaks, 
                f = exp(-breaks[1:(n-1)]), 
                df = -exp(-breaks[1:(n-1)]), 
                ddf.sup = exp(-breaks[1:(n-1)]))

# compare envelope to target function
curve(exp(-x), from = xmin, to = xmax)
curve(q2$e(x), from = xmin, to = xmax, add = TRUE, col = 4, n = 1e3)


#
# example 3: f(x) = sin(x)
#

# build envelope
q3 = bound.quad(breaks = breaks, 
                f = sin(breaks[1:(n-1)]) + 1, 
                df = cos(breaks[1:(n-1)]), 
                ddf.sup = sapply(1:(n-1), function(i) {
                  max(-sin(c(breaks[i:(i+1)])))
                }))

# build envelope with shifted anchors
anchors = breaks[1:(length(breaks)-1)] + diff(breaks)/2
q4 = bound.quad(breaks = breaks, 
                f = sin(anchors) + 1, 
                df = cos(anchors), 
                ddf.sup = sapply(1:(n-1), function(i) {
                  extrema = floor(2*breaks[i + 0:1]/pi) * pi / 2
                  extrema = extrema[
                    (extrema >= breaks[i]) & extrema <= breaks[i+1]
                  ]
                  max(-sin(c(extrema, breaks[i + 0:1])) )
                }), 
                anchors = anchors)

# compare envelope to target function
curve(sin(x) + 1, from = xmin, to = xmax)
curve(q3$e(x), from = xmin, to = xmax, add = TRUE, col = 4, n = 1e3)
curve(q4$e(x), from = xmin, to = xmax, add = TRUE, col = 2, n = 1e3)

