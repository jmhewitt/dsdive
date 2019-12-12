data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

spec = list(beta.sd = rep(1, 3), beta.absmax = 5, lambda.sd = 1, 
            sub.tx.mean = 0, sub.tx.sd = 1, surf.tx.mean = 0, surf.tx.sd = 1)

cl = snow::makeCluster(spec = 1, type = 'SOCK')

cfg = makeCompleteDistributed(dives = list(sim, sim), 
                              depth.bins = list(depth.bins, depth.bins), 
                              cl = cl, init = params, priors = spec, 
                              model = 'conditional')

x = dsdive.fit.gibbs.cfg(cfg = cfg, it = 2, verbose = TRUE, init = params, 
                         priors = spec, sigma = list(diag(3), diag(7)))

snow::stopCluster(cl)

detach(dive.sim$params)
detach(dive.sim)
