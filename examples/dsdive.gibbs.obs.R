data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

t.stages = sim$times[c(FALSE,diff(sim$stages)==1)]
  
tstep = diff(sim.obs$times[1:2])

obstx.mat = lapply(1:3, function(s) {
  dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta, 
                      lambda = lambda, s0 = s, tstep = tstep, 
                      include.raw = TRUE, delta = 1e-10)
})

lambda.priors = list(
  c(4, 2 + 2/3), c(4, 13 + 1/3), c(4, 5)
)

beta.priors = list(
  c(31.5, 3.5), c(3.5, 31.5)
)

T1.prior.params = c(25, .04)
T2.prior.params = c(56, .06)

fit = dsdive.gibbs.obs(
  dsobs.list = list(sim.obs), t.stages.list = list(t.stages), 
  beta.init = beta, lambda.init = lambda, verbose = TRUE, maxit = 2, 
  pi1.prior = beta.priors[[1]], pi2.prior = beta.priors[[2]], 
  lambda1.prior = lambda.priors[[1]], lambda2.prior = lambda.priors[[2]], 
  lambda3.prior = lambda.priors[[3]], tstep = tstep, depth.bins = depth.bins, 
  T1.prior.params = T1.prior.params, T2.prior.params = T2.prior.params, 
  max.width = 100, max.width.offset = 30, t0.prior.params = c(1,1), 
  tf.prior.params = c(1,1), offsets = 0, offsets.tf = 0, warmup = 1, 
  delta = 1e-10)


detach(dive.sim$params)
detach(dive.sim)
