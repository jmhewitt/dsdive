data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

# extract time between observations
tstep = diff(sim.obs$times[1:2])

# get true stage transition times
t.stages = sim$times[c(FALSE, diff(sim$stages)==1)]

# build probability matrices for observations
obstx.mat = lapply(1:3, function(s) {
  dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta, lambda = lambda, 
                      s0 = s, tstep = tstep, include.raw = TRUE)
})

# compute likelihood of observations, given model parameters
ld = dsdive.obsld(dsobs.list = list(sim.obs),
                  t.stages.list = list(t.stages), 
                  P.raw = obstx.mat)

detach(dive.sim$params)
detach(dive.sim)
