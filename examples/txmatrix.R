data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

m = dsdive.tx.matrix(t0 = 60, depth.bins = depth.bins, beta = beta, 
                     lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx)
