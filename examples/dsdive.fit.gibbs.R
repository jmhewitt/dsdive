data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

x = dsdive.fit.gibbs(depths = sim.obs$depths, times = sim.obs$times, 
                     durations = NULL, stages = NULL, depth.bins = depth.bins, 
                     t0.dive = sim.obs$times[1], it = 50, verbose = TRUE, 
                     init = params, priors.sd = rep(1,12), 
                     priors.mean = rep(0, 12))

detach(dive.sim$params)
detach(dive.sim)