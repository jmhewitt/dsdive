data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

# true stage transition times for the dive
t.stages = sim$times[c(FALSE, diff(sim$stages)==1)]

# uniformized transition rate
rate.unif = max(outer(lambda, 2*depth.bins[,2], '/'))

# time between observations
tstep = diff(sim.obs$times[1:2])

# probability transition matrix for observations
P.raw = lapply(1:3, function(s) {
  dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta, 
                      lambda = lambda, s0 = s, tstep = tstep, 
                      include.raw = TRUE)
})

# probability transition matrix for uniformized DTMC
P.tx = lapply(1:3, function(s) {
  dsdive.tx.matrix.uniformized(depth.bins = depth.bins, beta = beta, 
                               lambda = lambda, s0 = s, 
                               rate.uniformized = rate.unif)
})

# impute dive trajectory
x = dsdive.impute(
  depths = sim.obs$depths, times = sim.obs$times, t.stages = t.stages, 
  rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, n.bins = nrow(depth.bins), 
  max.tx = 100)

# overlay trajectory over observations
pl = plot(x = sim.obs, depth.bins = depth.bins, stages = sim.obs$stages, 
          errorbars = TRUE, imputed.list = x)

detach(dive.sim$params)
detach(dive.sim)
