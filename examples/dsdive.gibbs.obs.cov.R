data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)
library(parallel)
library(Rdsm)

cl = makeCluster(2, 'SOCK')
clusterEvalQ(cl, library(dsdive))
clusterEvalQ(cl, library(Rdsm))

mgrinit(cl)

t.stages = sim$times[c(FALSE,diff(sim$stages)==1)]
  
tstep = diff(sim.obs$times[1:2])

obstx.mat = lapply(1:3, function(s) {
  dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta, 
                      lambda = lambda, s0 = s, tstep = tstep, 
                      include.raw = TRUE, delta = 1e-10)
})

lambda.priors = list(
  list(mu = 0, sd = 3),
  list(mu = 0, sd = 3),
  list(mu = 0, sd = 3)
)

beta.priors = list(
  list(mu = 0, sd = 3),
  list(mu = 0, sd = 3)
)

T1.prior.params = c(25, .04)
T2.prior.params = c(56, .06)

dsobs.list = list(sim.obs, sim.obs)
t.stages.list = list(t.stages, t.stages)

beta.init = list(
  c(qlogis(.5), 0),
  c(0, 1)
)

alpha.init = c(beta.init, list(c(-1,-1)))

covs = data.frame(x1 = c(.5, 1), x2 = c(0, .3))

pi.formula = ~x1
lambda.formula = ~x1:x2

fit = dsdive.gibbs.obs.cov(
  dsobs.list = dsobs.list, t.stages.list = t.stages.list, 
  beta.init = beta.init, alpha.init = alpha.init, verbose = TRUE, maxit = 1, 
  beta1.prior = beta.priors[[1]], beta2.prior = beta.priors[[2]], 
  alpha1.prior = lambda.priors[[1]], alpha2.prior = lambda.priors[[2]], 
  alpha3.prior = lambda.priors[[3]], tstep = tstep, depth.bins = depth.bins, 
  T1.prior.params = T1.prior.params, T2.prior.params = T2.prior.params, 
  max.width = 100, max.width.offset = 30, t0.prior.params = c(1,1), 
  tf.prior.params = c(1,1), offsets = 0, offsets.tf = 0, warmup = 1, 
  covs = covs, pi.formula = pi.formula,  lambda.formula = lambda.formula, 
  cl = cl, optim.maxit = 1, delta = 1e-10)


detach(dive.sim$params)
detach(dive.sim)

stopCluster(cl)

