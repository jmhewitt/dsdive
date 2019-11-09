data('dive.sim')
attach(dive.sim)

out = fit.fixedimputation(par = params, imputed.list = list(imputed), 
                          it = 5, burn = 1, verbose = TRUE)
