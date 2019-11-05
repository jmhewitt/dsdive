data('dive.sim')
attach(dive.sim)
attach(params)

br = dsdive.bridgesample(depth.bins = depth.bins, d0 = sim.obs$depths[1], 
                         d0.last = NULL, df = sim.obs$depths[2], beta = beta, 
                         lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
                         t0 = sim.obs$times[1], tf = sim.obs$times[2], s0 = 1,
                         verbose = TRUE, precompute.bridges = TRUE)
