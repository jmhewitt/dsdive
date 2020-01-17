# build and save a simulated dive trajectory

library(dsdive)

#
# define domain and parameters
#

# load depth bins
depth.bins = read.csv('data_gen/depth_template.csv')

# define transition parameters
beta = c(.8, .2)
lambda = c(1, .3, .5)


#
# simulate and observe a dive
#

# simulate dive
x = dsdive.fwdsample.dive(depth.bins = depth.bins, beta = beta, lambda = lambda, 
                          t0 = 0, steps.max = 1e3, T1 = 15*60, T2 = 12*60)

# observe dive at regular time intervals
obs = dsdive.observe(depths = x$depths, times = x$times, stages = x$stages, 
                     t.obs = seq(from = 0, to = max(x$times) + 60, by = 60))

plot(x = x, depth.bins = depth.bins, dsobs = obs)

#
# package and export dive
#

dive.sim = list(
  params = list(beta = beta, lambda = lambda),
  depth.bins = depth.bins,
  sim = x,
  sim.obs = obs
)

save(dive.sim, file = 'data/dive.sim.RData')
