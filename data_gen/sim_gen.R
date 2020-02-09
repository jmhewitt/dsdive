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

# mean and sd. for time in stage
T1.mean = 10
T1.sd = 2
T2.mean = 15
T2.sd = 2

# convert mean and sd. for time in stage to gamma shape/rate parameters
T1.params = c(T1.mean^2/T1.sd^2, T1.mean/T1.sd^2)
T2.params = c(T2.mean^2/T2.sd^2, T2.mean/T2.sd^2)

# sample stage durations
T1 = 60 * rgamma(n = 1, shape = T1.params[1], rate = T1.params[2])
T2 = 60 * rgamma(n = 1, shape = T2.params[1], rate = T2.params[2])

# simulate dive
x = dsdive.fwdsample.dive(depth.bins = depth.bins, beta = beta, lambda = lambda, 
                          t0 = 0, steps.max = 1e3, T1 = T1, T2 = T2)

# observe dive at regular time intervals
obs = dsdive.observe(depths = x$depths, times = x$times, stages = x$stages, 
                     t.obs = seq(from = 0, to = max(x$times) + 60, by = 5*60))

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
