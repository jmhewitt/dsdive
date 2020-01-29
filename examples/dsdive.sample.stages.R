data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

# define prior parameters
T1.prior.params = c(15, 1/60)
T2.prior.params = c(15, 1/60)

# extract true stage transition times
t.stages.truth = sim$times[c(FALSE, diff(sim$stages)==1)]
t.stages = c(5*60, 50*60)

# sample stage transition times from full conditional posteriors
res = dsdive.sample.stages(depths = sim$depths, times = sim$times, 
                           t.stages = t.stages, beta = beta, lambda = lambda, 
                           depth.bins = depth.bins, 
                           T1.prior.params = T1.prior.params, 
                           T2.prior.params = T2.prior.params, max.width = 100, 
                           debug = TRUE)

rbind(original = t.stages, new = res$t.stages)

par(mar = c(5, 5, 4, 2) + .1)
# plot un-normalized full conditional posterior for T^(1) | T^(2), \cdot
curve(exp(res$debug$lp(x = x*60, stage.ind = 1, t.stages = t.stages)), 
      from = sim$times[1]/60, to = t.stages[2]/60, n = 1e2,
      xlab = expression(T^(1)), ylab = expression(f(T^(1)~'|'~T^(2), ...)))
# overlay sampling envelope
curve(res$debug$q1$e(x = x*60), from = sim$times[1]/60, to = t.stages[2]/60, 
      add = TRUE, col = 2, n = 1e2)
# overlay true stage 1->2 transition time
abline(v = t.stages.truth[1]/60, lty = 3)


par(mar = c(5, 5, 4, 2) + .1)
# plot un-normalized full conditional posterior for T^(2) | T^(1), \cdot
curve(exp(res$debug$lp(x = x*60, stage.ind = 2, 
                       t.stages = c(res$t.stages[1], t.stages[2]))), 
      from = t.stages[1]/60, to = sim$times[length(sim$times)]/60,
      xlab = expression(T^(2)), ylab = expression(f(T^(2)~'|'~T^(1), ...)))
# overlay sampling envelope
curve(res$debug$q2$e(x = x*60), 
      from = t.stages[1]/60, to = sim$times[length(sim$times)]/60,
      add = TRUE, col = 2, n = 1e2)
# overlay true stage 2->3 transition time
abline(v = t.stages.truth[2]/60, lty = 3)


#
# make comparison plots
#

dive.orig = dsdive.augment.trajectory(depths = sim$depths, times = sim$times, 
                                     t.stages = t.stages)
dive.orig = dsdive.simplify.trajectory(depths = dive.orig$depths, 
                                      times = dive.orig$times, 
                                      stages = dive.orig$stages)

dive.new = dsdive.augment.trajectory(depths = sim$depths, times = sim$times, 
                                     t.stages = res$t.stages)
dive.new = dsdive.simplify.trajectory(depths = dive.new$depths, 
                                      times = dive.new$times, 
                                      stages = dive.new$stages)

# examine original dive
pl1 = plot(x = dive.orig, depth.bins = depth.bins)
# compare with resampled stage transition times
pl2 = plot(x = dive.new, depth.bins = depth.bins)
# compare with truth
pl3 = plot(x = sim, depth.bins = depth.bins)

ggpubr::ggarrange(pl1, pl2, pl3, nrow = 3, ncol = 1, 
                  labels = c('Original:', 'Resampled:', 'Truth:'))

detach(dive.sim$params)
detach(dive.sim)
