data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

ld = dsdive.ld(depths = sim$depths, durations = sim$durations, 
               times = sim$times, stages = sim$stages, beta = beta, 
               lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
               depth.bins = depth.bins, t0.dive = 0, t.stage2 = NA, 
               model = 'conditional')

detach(dive.sim$params)
detach(dive.sim)
