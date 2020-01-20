data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

r.unif = max(outer(lambda, 2 * depth.bins$halfwidth, '/'))

m = dsdive.tx.matrix.uniformized(depth.bins = depth.bins, beta = beta, 
                                 lambda = lambda, s0 = 2, 
                                 rate.uniformized = r.unif)

N = 10

B = lapply(1:N, function(i) m)

x = ffbs.segment(B = B, x0 = 5, xN = 10, N = N)

detach(dive.sim$params)
detach(dive.sim)
