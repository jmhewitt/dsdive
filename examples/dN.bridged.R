data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

r.unif = max(outer(lambda, 2 * depth.bins$halfwidth, '/'))

m = dsdive.tx.matrix.uniformized(depth.bins = depth.bins, beta = beta, 
                                 lambda = lambda, s0 = 2, 
                                 rate.uniformized = r.unif)

N.max = 20

dN = dN.bridged(B = m, x0 = 3, xN = 5, N.max = N.max, rate.uniformized = r.unif, 
                t = 60, log = FALSE)

plot(0:N.max, dN)

detach(dive.sim$params)
detach(dive.sim)
