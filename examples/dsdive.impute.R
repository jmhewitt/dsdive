data('dive.sim')
attach(dive.sim)
attach(params)

imputed = dsdive.impute(depth.bins = depth.bins, depths = sim.obs$depths, 
                        times = sim.obs$times, s0 = 1, beta = beta, 
                        lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx,
                        t0.dive = 0)

detach(params)
detach(dive.sim)