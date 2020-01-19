data('dive.sim')
attach(dive.sim)
attach(params)

library(snow)
library(parallel)

stopCluster(cl)
cl = makeCluster(detectCores() - 1, 'SOCK')

devtools::document()
res = dsdive.fit.tpois(dives.obs = list(list(dive = sim.obs, 
                                             depth.bins = depth.bins)), 
                 cl = cl, 
                 beta.init = c(.5,.5), lambda.init = lambda, verbose = TRUE, 
                 T1.prior = function(x) dgamma(x, 36, 3, log = TRUE), 
                 T2.prior = function(x) dgamma(x, 36, 2, log = TRUE), 
                 pi1.prior = c(1,1), pi2.prior = c(1,1), lambda1.prior = c(2,1), 
                 lambda2.prior = c(2,1), lambda3.prior = c(2,1), it = 5, 
                 inflation.factor.lambda = 1.1)

stopCluster(cl)

res$theta

detach(params)
detach(dive.sim)