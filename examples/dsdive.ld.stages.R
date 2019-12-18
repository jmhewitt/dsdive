data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

dens = dsdive.ld.stages(breaks = c(25, 50), fixed.ind = 2, beta = beta, 
                        lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
                        depths = sim$depths, durations = sim$durations, 
                        times = sim$times, depth.bins = depth.bins, t0.dive = 0, 
                        t.stage2 = sim$times[min(which(sim$stages==2))],
                        model = 'conditional')

detach(dive.sim$params)
detach(dive.sim)
