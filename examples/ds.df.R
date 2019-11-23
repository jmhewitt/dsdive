data('dive.sim')
attach(dive.sim)

df = ds.df(depths = sim.obs$depths, times = sim.obs$times, 
           depth.bins = depth.bins, stages = NULL)

detach(dive.sim)
