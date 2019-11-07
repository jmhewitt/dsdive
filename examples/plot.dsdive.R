data('dive.sim')
attach(dive.sim)

# basic plot
pl = plot(x = sim, depth.bins = depth.bins)

# plot with observations overlaid
pl2 = plot(x = sim, depth.bins = depth.bins, dsobs = sim.obs)

# plot with imputed trajectory overlaid
pl3 = plot(x = sim, depth.bins = depth.bins, dsobs = sim.obs, 
           imputed.list = imputed)
