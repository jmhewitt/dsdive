# define range for the envelope
xmin = 1
xmax = 10

# define envelope segments
breaks = seq(from = xmin, to = xmax, length.out = 10)

# define midpoints of segments
anchors = breaks[1:(length(breaks)-1)] + diff(breaks)/2

# get number of segments
n = length(breaks)

#
# example 1: f(x) = exp(-x); log f(x) = -x
#

e.approx = envelope.approx(breaks = breaks, anchors = anchors, 
                           lp.breaks = -breaks, lp.anchors = -anchors)
