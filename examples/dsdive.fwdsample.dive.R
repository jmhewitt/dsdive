data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

# simulate dive
x = dsdive.fwdsample.dive(depth.bins = depth.bins, beta = beta, 
                          lambda = lambda, t0 = 0, steps.max = 1e3, 
                          T1 = 15*60, T2 = 12*60)

pl = plot(x = x, depth.bins = depth.bins)

detach(dive.sim$params)
detach(dive.sim)
