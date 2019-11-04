data('dive.sim')
attach(dive.sim)

obs = dsdive.observe(depths = sim$depths, times = sim$times, 
                     stages = sim$stages, 
                     t.obs = seq(from = 0, to = max(sim$times)+60, by = 60))
