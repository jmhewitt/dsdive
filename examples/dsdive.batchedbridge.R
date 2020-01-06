data('dive.sim')
attach(dive.sim)
attach(params)

depth.bins = depth.bins[1:16,]

beta = c(.8,.2)

br = dsdive.batchedbridge(depth.bins = depth.bins, depths = c(1,5, 5, 10), 
                          d0.last = NULL, times = c(0,60, 62,120), tx.max = 100,
                          beta = beta, lambda = lambda, sub.tx = sub.tx, 
                          surf.tx = surf.tx, s0 = 1, sf = 1, 
                          inflation.factor.lambda = 1.1, verbose = FALSE, 
                          lambda.max = NULL, t0.dive = 0, 
                          t.stages = c(500, 1e3), model = 'logit')

detach(params)
detach(dive.sim)
