data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

o = dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta, lambda = lambda, 
                        s0 = 1, tstep = 300, include.raw = TRUE)


detach(dive.sim$params)
detach(dive.sim)
