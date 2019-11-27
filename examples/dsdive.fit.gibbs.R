data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

spec = list(beta.sd = rep(1, 3), beta.absmax = 5, lambda.sd = 1, 
            sub.tx.mean = 0, sub.tx.sd = 1, surf.tx.mean = 0, surf.tx.sd = 1)

x = dsdive.fit.gibbs(depths = sim.obs$depths, times = sim.obs$times, 
                     durations = NULL, stages = NULL, depth.bins = depth.bins, 
                     t0.dive = sim.obs$times[1], it = 50, verbose = TRUE, 
                     init = params, priors = spec)

detach(dive.sim$params)
detach(dive.sim)