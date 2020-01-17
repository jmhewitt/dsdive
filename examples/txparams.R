data('dive.sim')
attach(dive.sim)
attach(params)

dsdive.tx.params(depth.bins = depth.bins, d0 = 3, s0 = 1, 
                 beta = beta, lambda = lambda)

detach(params)
detach(dive.sim)
