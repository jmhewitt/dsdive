data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

T1.logprior = function(x) dgamma(x, shape = 15, rate = 1, log = TRUE)
T2.logprior = function(x) dgamma(x, shape = 15, rate = 1, log = TRUE)

t.stages = c(5, 30)*60

augmented = dsdive.augment.trajectory(depths = sim$depths, times = sim$times, 
                                      t.stages = t.stages)

dive.new = dsdive.sample.stages(depths = sim$depths, times = sim$times, 
                                t.stages = t.stages, beta = beta, 
                                lambda = lambda, depth.bins = depth.bins, 
                                T1.logprior = T1.logprior, 
                                T2.logprior = T2.logprior, 
                                t1.sd = 1, t2.sd = 1)

rbind(t.stages, dive.new$t.stages)

pl1 = plot(x = augmented, depth.bins = depth.bins)
pl2 = plot(x = dive.new$augmented, depth.bins = depth.bins)

ggpubr::ggarrange(pl1, pl2, nrow = 2, ncol = 1)

detach(dive.sim$params)
detach(dive.sim)
