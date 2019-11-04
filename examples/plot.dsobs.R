data('dive.sim')
attach(dive.sim)

# basic plot
pl = plot(x = sim.obs, depth.bins = depth.bins)

# plot with stages overlaid
pl2 = plot(x = sim.obs, depth.bins = depth.bins, stages = sim.obs$stages)
