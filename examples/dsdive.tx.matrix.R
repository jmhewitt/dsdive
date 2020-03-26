data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

m = dsdive.tx.matrix(t0 = 60, depth.bins = depth.bins, beta = beta, 
                     lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx,
                     t0.dive = 0, t.stage2 = NA, model = 'conditional')

detach(dive.sim$params)
detach(dive.sim)