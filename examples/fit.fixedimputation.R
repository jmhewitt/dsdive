data('dive.sim')
attach(dive.sim)

out = fit.fixedimputation(par = params, imputed.list = list(imputed), 
                          it = 1, burn = 1, verbose = TRUE, t0.dive = 0)

detach(dive.sim)