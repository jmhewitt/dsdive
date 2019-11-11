data('dive.sim')
attach(dive.sim)

x = dsdive.stochit(depths = sim.obs$depths, t = sim.obs$times, 
                   depth.bins = depth.bins, par = params, 
                   priors = rep(5,12), max.resample = 5, max.it = 1e2, tol = .1, 
                   n.paths = 1e1, n.par = 1e1, verbose = TRUE, 
                   control.mh = list(burn_in = 1e2, quiet = FALSE))

detach(dive.sim)