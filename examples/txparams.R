data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

dsdive.tx.params(t0 = 0, depth.bins = depth.bins, d0 = 0, 
                 d0.last = NULL, s0 = 1, beta = beta, lambda = lambda, 
                 sub.tx = sub.tx, surf.tx = surf.tx, t0.dive = 0, t.stage2 = NA)

dsdive.tx.params(t0 = 0, depth.bins = depth.bins, d0 = 60, 
                 d0.last = 61, s0 = 1, beta = beta, lambda = lambda, 
                 sub.tx = sub.tx, surf.tx = surf.tx, t0.dive = 0, t.stage2 = NA)

dsdive.tx.params(t0 = 0, depth.bins = depth.bins, d0 = 60, 
                 d0.last = 61, s0 = 2, beta = beta, lambda = lambda, 
                 sub.tx = sub.tx, surf.tx = surf.tx, t0.dive = 0, t.stage2 = NA)

dsdive.tx.params(t0 = 0, depth.bins = depth.bins, d0 = 60, 
                 d0.last = 61, s0 = 2, beta = beta, lambda = lambda, 
                 sub.tx = sub.tx, surf.tx = surf.tx, t0.dive = 0, 
                 shift.params = c(50, 1.125, 3, .5), t.stage2 = NA)

detach(dive.sim$params)
detach(dive.sim)