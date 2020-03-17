data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

# extract true stage transition times
t.stages = sim$times[c(FALSE, diff(sim$stages)==1)]

# build probability matrices for observations
tstep = diff(sim.obs$times[1:2])
P.raw = lapply(1:3, function(s) {
  dsdive.obstx.matrix(depth.bins = depth.bins, beta = params$beta, 
                      lambda = params$lambda, s0 = s, tstep = tstep, 
                      include.raw = TRUE)
})

# sample new offsets
res = dsdive.obs.sample.offsets(
  dsobs.aligned = sim.obs, dsobs.unaligned = sim.obs, offset = 0, 
  t.stages = t.stages, P.raw = P.raw, depth.bins = depth.bins, tstep = tstep,
  prior.params = c(10,10), max.width = 30, debug = TRUE, sample.start = TRUE, 
  offset.tf = 0)


par(mar = c(5, 5, 4, 2) + .1)
xmin = -tstep
xmax = tstep
# plot un-normalized full conditional posterior for offset
curve(exp(res$debug$lp(eps = x)), 
      from = xmin, to = xmax, n = 1e3, col = 'grey60',
      xlab = expression(epsilon), ylab = expression(f(epsilon~'|'~...)))
# overlay sampling envelope
curve(exp(res$debug$q1$e.log(x = x)), from = xmin, to = xmax, 
      add = TRUE, col = 2, n = 1e3)

detach(dive.sim$params)
detach(dive.sim)
