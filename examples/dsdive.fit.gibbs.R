data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

# simulate dive
devtools::document()
x = dsdive.fit.gibbs(depths = sim.obs$depths, times = sim.obs$times, 
                     durations = NULL, stages = NULL, depth.bins = depth.bins, 
                     t0.dive = sim.obs$times[1], it = 50, verbose = TRUE, 
                     init = params, sigma = diag(12))

plot(x = sim, depth.bins = depth.bins, dsobs = sim.obs, 
     imputed.list = x$trace.imputed, imputed.alpha = .03)

detach(dive.sim$params)
detach(dive.sim)