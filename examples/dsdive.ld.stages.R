data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

dens = dsdive.ld.stages(breaks = c(200, 400), fixed.ind = 2, beta = beta, 
                        lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
                        depths = sim$depths, durations = sim$durations, 
                        times = sim$times, depth.bins = depth.bins)