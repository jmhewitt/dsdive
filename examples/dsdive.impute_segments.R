data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

t.sbreaks = sim$times[c(FALSE, diff(sim$stages)==1)]

x = dsdive.impute_segments(depth.bins = depth.bins, 
                           depths = sim.obs$depths, 
                           times = sim.obs$times, beta = beta, 
                           lambda = lambda, s0 = 1, 
                           inflation.factor.lambda = 1.1, 
                           method.N = 'truncpois', N.max = 100, 
                           t.sbreaks = t.sbreaks )

pl = plot(x = sim.obs, depth.bins = depth.bins, stages = sim.obs$stages, 
          errorbars = TRUE, imputed.list = x, imputed.alpha = .6)

detach(dive.sim$params)
detach(dive.sim)
