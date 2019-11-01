data('dive.sim')
attach(dive.sim)

stages.est = c(rep(1,4), rep(2,15), rep(3, 7))

est = dsdive.heurest(depths = sim.obs$depths, times = sim.obs$times, 
                     stages.est = stages.est, depth.bins = depth.bins, 
                     sub.tx1 = params$sub.tx[1])
