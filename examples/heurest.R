data('dive.sim')
attach(dive.sim)

spec = list(beta.sd = rep(1, 3), beta.absmax = 5, lambda.sd = 1, 
            sub.tx.mean = 0, sub.tx.sd = 1, surf.tx.mean = 0, surf.tx.sd = 1,
            sub.tx = params$sub.tx, surf.tx = params$surf.tx)

est = dsdive.heurest(depths = sim.obs$depths, times = sim.obs$times, 
                     stages.est = sim.obs$stages, depth.bins = depth.bins,
                     t0.dive = 0, priors = spec, model = 'conditional')

detach(dive.sim)