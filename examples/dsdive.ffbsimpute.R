data('dive.sim')
attach(dive.sim)
attach(params)



t.stages = seq(from = min(sim.obs$times), to = max(sim.obs$times), 
               length.out = 4)[2:3]

beta = c(.8,.2)
devtools::document()

imputed = dsdive.ffbsimpute(depth.bins = depth.bins, depths = sim.obs$depths, 
                            times = sim.obs$times, s0 = 1, beta = beta, 
                            lambda = lambda, sub.tx = sub.tx, 
                            surf.tx = surf.tx, 
                            inflation.factor.lambda = 1.1, t0.dive = 0, 
                            t.stages = t.stages, model = 'logit', 
                            times.cond = NULL)

all(imputed$times[imputed$stages==1] < t.stages[1])
all(imputed$times[imputed$stages==2] < t.stages[2])


plot(x = sim.obs, depth.bins = depth.bins, errorbars = TRUE, 
     imputed.list = sim, imputed.alpha = .8)

detach(params)
detach(dive.sim)
