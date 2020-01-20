data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

t.stages = sim$times[diff(sim$stages) == 1]

T1.prior = function(x) dgamma(x, shape = 15, rate = 1, log = TRUE)
T2.prior = function(x) dgamma(x, shape = 15, rate = 1, log = TRUE)

dsdive.sample.stages(depths = sim$depths, durations = sim$durations, 
                     times = sim$times, t.stages = t.stages, beta = beta, 
                     lambda = lambda, depth.bins = depth.bins, 
                     T1.prior = T1.prior, T2.prior = T2.prior)

detach(dive.sim$params)
detach(dive.sim)
