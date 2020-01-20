data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

t.stages = sim$times[c(FALSE, diff(sim$stages)==1)]

x = dsdive.impute(depth.bins = depth.bins, depths = sim.obs$depths, 
                  times = sim.obs$times, beta = beta, lambda = lambda, 
                  inflation.factor.lambda = 1.1, t.stages = t.stages, 
                  method.N = 'truncpois', N.max = 100)

x$durations = c(x$durations, Inf)

pl = plot(x = sim.obs, depth.bins = depth.bins, stages = sim.obs$stages, 
          errorbars = TRUE, imputed.list = x, imputed.alpha = .6)

pl

detach(dive.sim$params)
detach(dive.sim)