# build and save a simulated dive trajectory

library(dsdive)

#
# define domain and parameters
#

# define depth domain
bin.width = 5
bin.centers = seq(from = bin.width/2, to = 1e3, by = bin.width)
depths = cbind(center = bin.centers, halfwidth = bin.width/2)

# define transition parameters
beta = matrix(c(2.5,   0, -2.5, 
                0, -.5, 0), 
              nrow = 2, byrow = TRUE)
lambda = c(.5,1,.5)
sub.tx = c(10, 15)
surf.tx = sub.tx



#
# simulate and observe a dive
#

# simulate dive
x = dsdive.fwdsample(depth.bins = depths, d0 = 1, beta = beta, 
                     lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
                     t0 = 0, tf = Inf, steps.max = 1e5, dur0 = NULL, 
                     nsteps = NULL, s0 = 1, t0.dive = 0, t.stage2 = NA)

# observe dive at regular time intervals
obs = dsdive.observe(depths = x$depths, times = x$times, stages = x$stages,
                     t.obs = seq(from = 0, to = max(x$times)+60, by = 1*60))

plot(x = x, depth.bins = depths, dsobs = obs)

# impute trajectory
imputed = dsdive.fastimpute(M = 1, depth.bins = depths, depths = obs$depths, 
                            times = obs$times, s0 = 1, beta = beta, 
                            lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
                            t0.dive = 0)

plot(x = x, depth.bins = depths, dsobs = obs, imputed.list = imputed, 
     imputed.alpha = .3)


#
# package and export dive
#

dive.sim = list(
  params = list(beta = beta, lambda = lambda, sub.tx = sub.tx, 
                surf.tx = surf.tx),
  depth.bins = depths,
  sim = x,
  sim.obs = obs,
  imputed = imputed
)

save(dive.sim, file = 'data/dive.sim.RData')
