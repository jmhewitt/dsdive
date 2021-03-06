data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

# define prior parameters
T1.prior.params = c(15, 1/60)
T2.prior.params = c(15, 1/60)

# extract true stage transition times
s12.inds = which(sim$stages<3)
t.stages.truth = sim$times[c(FALSE, diff(sim$stages)==1)]
t.stages = c(20,40)*60#t.stages.truth

# extract times at which we observe depth bin transitions
depth.bin.tx = sim.obs$times[c(FALSE, diff(sim.obs$depths) != 0)]

# build probability matrices for observations
tstep = diff(sim.obs$times[1:2])
P.raw = lapply(1:3, function(s) {
  dsdive.obstx.matrix(depth.bins = depth.bins, beta = params$beta, 
                      lambda = params$lambda, s0 = s, tstep = tstep, 
                      include.raw = TRUE, delta = 1e-10)
})

# sample stage transition times from full conditional posteriors
res = dsdive.obs.sample.stages(
  depths = sim.obs$depths, times = sim.obs$times, t.stages = t.stages, 
  P.raw = P.raw, T.range = range(sim$times), depth.bins = depth.bins, 
  T1.prior.params = T1.prior.params, T2.prior.params = T2.prior.params, 
  max.width = 100, debug = TRUE)
  
rbind(original = t.stages, new = res$t.stages)

par(mar = c(5, 5, 4, 2) + .1)
xmin = sim$times[1]/60
xmax = t.stages[2]/60
# plot un-normalized full conditional posterior for T^(1) | T^(2), \cdot
curve(res$debug$lp(x = x*60, stage.ind = 1, t.stages = t.stages), 
      from = xmin, to = xmax, n = 1e3, col = 'grey60',
      xlab = expression(T^(1)), ylab = expression(ln~f(T^(1)~'|'~T^(2), ...)))
# overlay observed depth bin transition times
abline(v = depth.bin.tx/60, lty = 3)
# overlay sampling envelope
curve(res$debug$q1$e.log(x = x*60), from = xmin, to = xmax, 
      add = TRUE, col = 2, n = 1e3)
# overlay true stage 1->2 transition time
abline(v = t.stages.truth[1]/60, lty = 2)


par(mar = c(5, 5, 4, 2) + .1)
xmin = res$t.stages[1]/60
xmax = sim$times[length(sim$times)]/60
# plot un-normalized full conditional posterior for T^(2) | T^(1), \cdot
curve(exp(res$debug$lp(x = x*60, stage.ind = 2, 
                       t.stages = c(res$t.stages[1], t.stages[2]))), 
      from = xmin, to = xmax, n = 1e3,
      xlab = expression(T^(2)), ylab = expression(f(T^(2)~'|'~T^(1), ...)))
# overlay sampling envelope
curve(res$debug$q2$e(x = x*60), from = xmin, to = xmax, 
      add = TRUE, col = 2, n = 1e3)
# overlay true stage 2->3 transition time
abline(v = t.stages.truth[2]/60, lty = 3)

detach(dive.sim$params)
detach(dive.sim)
