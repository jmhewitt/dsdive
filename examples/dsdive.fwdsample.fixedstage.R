data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

# simulate dive
x = dsdive.fwdsample.fixedstage(depth.bins = depth.bins, d0 = 1, 
                                beta = c(.8, .2), lambda = rep(1,3), t0 = 0, 
                                tf = 100, steps.max = 1e3, s0 = 1)

pl = plot(x = x, depth.bins = depth.bins)

detach(dive.sim$params)
detach(dive.sim)
