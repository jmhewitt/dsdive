data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

ld = dsdive.ldabc(beta = beta, lambda = lambda, sub.tx = sub.tx, 
                  surf.tx = surf.tx, depths.labels = depth.bins, 
                  steps.max = 1e3, N = 1e1, depths = sim.obs$depths, 
                  t = sim.obs$times, tries.max = 1e5, dump.state = TRUE, 
                  verbose = TRUE, n.samples = 5, eps = 1)
