data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

s1.inds = sim.obs$stages==1

x = dsdive.impute_segments(depth.bins = depth.bins, 
                           depths = sim.obs$depths[s1.inds], 
                           times = sim.obs$times[s1.inds], beta = beta, 
                           lambda = lambda, s0 = 1, 
                           inflation.factor.lambda = 1.1, method.N = 'exact', 
                           N.max = 100)

pl = plot(x = sim.obs, depth.bins = depth.bins, stages = sim.obs$stages, 
          errorbars = TRUE, imputed.list = x, imputed.alpha = .6)

detach(dive.sim$params)
detach(dive.sim)
