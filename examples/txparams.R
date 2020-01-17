data('dive.sim')
attach(dive.sim)

dsdive.tx.params(depth.bins = depth.bins, d0 = 3, s0 = 1, 
                 beta = c(.8, .2), lambda = rep(1,3))

detach(dive.sim)
