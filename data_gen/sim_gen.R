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
beta = matrix(c(2.5,   0, -1.5, 
                -.5, -.5, -.75), 
              nrow = 2, byrow = TRUE)
lambda = 1/c(3, 3, 3)
sub.tx = c(50, .02)
surf.tx = c(-10, 5e-3)


#
# simulate and observe a dive
#

# simulate dive
x = dsdive.fwdsample(depth.bins = depths, d0 = 1, beta = beta, 
                     lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
                     t0 = 0, tf = Inf, steps.max = 1e5, dur0 = NULL, 
                     nsteps = NULL, s0 = 1, t0.dive = 0)

# observe dive at regular time intervals
obs = dsdive.observe(depths = x$depths, times = x$times, stages = x$stages,
                     t.obs = seq(from = 0, to = max(x$times)+60, by = 1*60))

# impute trajectory
imputed = dsdive.impute(depth.bins = depths, depths = obs$depths, 
                        times = obs$times, s0 = 1, beta = beta, lambda = lambda, 
                        sub.tx = sub.tx, surf.tx = surf.tx, t0.dive = 0)


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
