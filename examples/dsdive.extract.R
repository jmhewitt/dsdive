data('dive.sim')
attach(dive.sim)
attach(params)

dsdive.extract(depths = sim$depths, times = sim$times, stages = sim$stages,
               durations = sim$durations, t0 = 0, tf = 60)

detach(params)
detach(dive.sim)