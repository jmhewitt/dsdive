data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

spec = list(beta.sd = rep(1, 3), beta.absmax = 5, lambda.sd = 1, 
            sub.tx.mean = 0, sub.tx.sd = 1, surf.tx.mean = 0, surf.tx.sd = 1)

cfg = makeImputedSingle(depth.bins = depth.bins, it = 5, 
                        depths = sim.obs$depths, times = sim.obs$times, 
                        init = params, priors = spec, 
                        inflation.factor.lambda = 1.1, t0.dive = 0, 
                        verbose = TRUE)

detach(dive.sim$params)
detach(dive.sim)
