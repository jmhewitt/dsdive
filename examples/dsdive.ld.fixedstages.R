data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

dsdive.ld.fixedstages(depths = sim$depths, durations = sim$durations, 
                      times = sim$times, stages = sim$stages, beta = beta, 
                      lambda = lambda, depth.bins = depth.bins)

detach(dive.sim$params)
detach(dive.sim)
