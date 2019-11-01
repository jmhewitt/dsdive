data('dive.sim')
attach(dive.sim)

# basic plot
pl = plot(x = sim.obs, depth.bins = depth.bins)

# plot with stages overlaid
stages = c(rep(1, 5), rep(2,15), rep(3,6))
pl2 = plot(x = sim.obs, depth.bins = depth.bins, stages = stages)
