data('dive.sim')
attach(dive.sim)

imputed = crawl.impute(depth.bins = depth.bins, depths = sim.obs$depths, 
                       times = sim.obs$times, N = 1, depths.impute = 'midpoint', 
                       tstep = 5)

detach(dive.sim)