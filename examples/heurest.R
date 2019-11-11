data('dive.sim')
attach(dive.sim)

est = dsdive.heurest(depths = sim.obs$depths, times = sim.obs$times, 
                     stages.est = sim.obs$stages, depth.bins = depth.bins, 
                     sub.tx1 = params$sub.tx[1], t0.dive = 0)

detach(dive.sim)