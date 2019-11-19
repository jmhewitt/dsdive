data('dive.sim')
attach(dive.sim)
attach(params)

br = dsdive.fastbridge(M = 1e2, depth.bins = depth.bins, d0 = sim.obs$depths[1], 
                       d0.last = NULL, df = sim.obs$depths[2], beta = beta, 
                       lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
                       t0 = sim.obs$times[1], tf = sim.obs$times[2], s0 = 1,
                       verbose = TRUE, t0.dive = 0, t.stage2 = NA)

detach(params)
detach(dive.sim)
