data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

spec = list(beta.sd = rep(1, 3), beta.absmax = 5, lambda.sd = 1, 
            sub.tx.mean = 0, sub.tx.sd = 1, surf.tx.mean = 0, surf.tx.sd = 1)

cl = snow::makeCluster(spec = 1, type = 'SOCK')

cfg = makeImputedDistributed(dives = list(sim.obs, sim.obs), 
                             depth.bins = list(depth.bins, depth.bins), 
                             cl = cl, init = params, priors = spec, it = 5, 
                             inflation.factor.lambda = 1.1)

snow::stopCluster(cl)

detach(dive.sim$params)
detach(dive.sim)
