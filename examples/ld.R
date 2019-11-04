data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

ld = dsdive.ld(depths = sim$depths, durations = sim$durations, 
               times = sim$times, stages = sim$stages, beta = beta, 
               lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
               depths.labels = depth.bins)
