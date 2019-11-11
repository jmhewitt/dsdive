data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

ld = dsdive.ldabc(beta = beta, lambda = lambda, sub.tx = sub.tx, 
                  surf.tx = surf.tx, depth.bins = depth.bins, 
                  steps.max = 1e3, N = 1e1, depths = sim.obs$depths, 
                  t = sim.obs$times, tries.max = 1e2, dump.state = TRUE, 
                  verbose = TRUE, n.samples = 5, eps = 1, t0.dive = 0)

detach(dive.sim$params)
detach(dive.sim)