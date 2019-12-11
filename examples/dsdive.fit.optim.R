data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

spec = list(beta.sd = rep(1, 3), beta.absmax = 5, lambda.sd = 1, 
            sub.tx.mean = 0, sub.tx.sd = 1, surf.tx.mean = 0, surf.tx.sd = 1)

cfg = makeCompleteSingle(depth.bins = depth.bins, durations = sim$durations, 
                         depths = sim$depths, times = sim$times, 
                         stages = sim$stages, init = params, priors = spec, 
                         t0.dive = 0)

x = dsdive.fit.optim(cfg = cfg, init = params, priors = spec, hessian = FALSE,
                     maxit = 1, method = 'Nelder-Mead')

detach(dive.sim$params)
detach(dive.sim)
