data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

dsdive.tx.params(t0 = 0, num.depths = length(depth.bins)-1, d0 = 0, 
                 d0.last = NULL, s0 = 1, beta = beta, lambda = lambda, 
                 sub.tx = sub.tx, surf.tx = surf.tx)

dsdive.tx.params(t0 = 0, num.depths = length(depth.bins)-1, d0 = 60, 
                 d0.last = 61, s0 = 1, beta = beta, lambda = lambda, 
                 sub.tx = sub.tx, surf.tx = surf.tx)

dsdive.tx.params(t0 = 0, num.depths = length(depth.bins)-1, d0 = 60, 
                 d0.last = 61, s0 = 2, beta = beta, lambda = lambda, 
                 sub.tx = sub.tx, surf.tx = surf.tx)