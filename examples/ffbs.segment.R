data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

r.unif = max(outer(lambda, 2 * depth.bins$halfwidth, '/'))

m = dsdive.tx.matrix.uniformized(depth.bins = depth.bins, beta = beta, 
                                 lambda = lambda, s0 = 2, 
                                 rate.uniformized = r.unif)

x = ffbs.segment(B = m$m, x0 = 5, xN = 10, N = 10)

detach(dive.sim$params)
detach(dive.sim)
