data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

t.stages = sim$times[c(FALSE,diff(sim$stages)==1)]
  
tstep = diff(sim.obs$times[1:2])

obstx.mat = lapply(1:3, function(s) {
  dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta, 
                      lambda = lambda, s0 = s, tstep = tstep, 
                      include.raw = TRUE)
})

lambda.priors = list(
  c(4, 2 + 2/3), c(4, 13 + 1/3), c(4, 5)
)

beta.priors = list(
  c(31.5, 3.5), c(3.5, 31.5)
)

fit = dsdive.gibbs.obs(
  dsobs.list = list(sim.obs), t.stages.list = list(t.stages), 
  beta.init = beta, lambda.init = lambda, verbose = TRUE, maxit = 1e3, 
  pi1.prior = beta.priors[[1]], pi2.prior = beta.priors[[2]], 
  lambda1.prior = lambda.priors[[1]], lambda2.prior = lambda.priors[[2]], 
  lambda3.prior = lambda.priors[[3]], tstep = tstep, depth.bins = depth.bins)


detach(dive.sim$params)
detach(dive.sim)
