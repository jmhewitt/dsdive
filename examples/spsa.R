data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

o = spsa(par = c(beta, log(lambda), log(sub.tx[-1]), surf.tx), 
          fn = function(theta) {
            r = dsdive.ldabc(beta = matrix(theta[1:6], nrow = 2), 
                             lambda = exp(theta[7:9]), 
                             sub.tx = c(sub.tx[1], exp(theta[10])), 
                             surf.tx = theta[11:12], 
                             depth.bins = depth.bins, 
                             steps.max = 1e3, N = 1e1, depths = sim.obs$depths, 
                             t = sim.obs$times, tries.max = 1e5, 
                             dump.state = TRUE, verbose = FALSE, n.samples = 1, 
                             eps = 1)$ld
            -r
          }, c = 2^(-6), a = 2^(-12), alpha = 1/3)

detach(dive.sim$params)
detach(dive.sim)