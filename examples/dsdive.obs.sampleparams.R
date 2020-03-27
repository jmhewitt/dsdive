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

x = dsdive.obs.sampleparams(dsobs.list = list(sim.obs), 
                            t.stages.list = list(t.stages), P.raw = obstx.mat, 
                            s0 = 1, depth.bins = depth.bins, beta = beta, 
                            lambda = lambda, lambda.priors.list = lambda.priors, 
                            beta.priors.list = beta.priors, 
                            tstep = diff(sim.obs$times[1:2]))

detach(dive.sim$params)
detach(dive.sim)
