data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

# get true stage transition times
t.stages = sim$times[c(FALSE, diff(sim$stages)==1)]

dsdive.suffstat(depths = sim$depths, times = sim$times, t.stages = t.stages, 
                depth.bins = depth.bins)

detach(dive.sim$params)
detach(dive.sim)
