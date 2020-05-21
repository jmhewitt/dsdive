data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

library(parallel)


#
# build interpolating functions
#

cl = makeCluster(spec = 3, type = 'SOCK')
clusterEvalQ(cl, library(dsdive))

interpolators = lapply(1:3, function(s0) {
  
  # support for directional preferences, speeds, and timesteps
  if(s0 == 1) {
    beta.seq = seq(.5, 1, length.out = 7)
    beta.seq[length(beta.seq)] = .999
    lambda.seq = seq(.1, 2, length.out = 7)
    tstep.seq = seq(0, 300, by = 100)
  } else if(s0 == 2) {
    beta.seq = .5
    lambda.seq = seq(.1, 2, length.out = 7)
    tstep.seq = seq(0, 300, by = 100)
  } else if(s0 == 3) {
    beta.seq = seq(0, .5, length.out = 7)
    beta.seq[1] = 1-.999
    lambda.seq = seq(.1, 2, length.out = 7)
    tstep.seq = seq(0, 300, by = 100)
  }
  
  interpolators = dsdive.obstx.matrix_interpolator(
    depth.bins = depth.bins, beta.seq = beta.seq, lambda.seq = lambda.seq,
    s0 = s0, tstep.seq = tstep.seq, m = 3, verbose = TRUE, cl = cl)
  
})

stopCluster(cl = cl)

# get true stage transition times
t.stages = sim$times[c(FALSE, diff(sim$stages)==1)]

# compute likelihood of observations, given model parameters
ld = dsdive.obsld_approx(dsobs.list = list(sim.obs,sim.obs),
                         t.stages.list = list(t.stages,t.stages), 
                         s0 = 1, sf = 3, beta = params$beta, 
                         lambda = params$lambda, P.interpolator = interpolators)

detach(dive.sim$params)
detach(dive.sim)
