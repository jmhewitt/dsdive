data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

# simulate dive
x = dsdive.fwdsample(depth.bins = depth.bins, d0 = 1, beta = beta, 
                     lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
                     t0 = 0, tf = Inf, steps.max = 1e5, dur0 = NULL, 
                     nsteps = NULL, s0 = 1)
