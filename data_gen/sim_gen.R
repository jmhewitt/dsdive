# build and save a simulated dive trajectory

library(dsdive)

#
# define domain and parameters
#

# define depth domain
max.depths = seq(from = 0, to = 1e3, by = 5)
num.depths = length(max.depths)
depth.bins = paste('[', max.depths[1:(num.depths-1)], ', ', max.depths[-1],
                   ')', sep = '')

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
x = dsdive.fwdsample(depths.labels = depth.bins, d0 = 0, beta = beta, 
                     lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
                     t0 = 0, tf = Inf, steps.max = 1e5, dur0 = NULL, 
                     nsteps = NULL, s0 = 1)

# observe dive at regular time intervals
obs = dsdive.observe(depths = x$depths, times = x$times, 
                     t.obs = seq(from = 0, to = max(x$times)+60, by = 1*60))


#
# package and export dive
#

dive.sim = list(
  params = list(beta = beta, lambda = lambda, sub.tx = sub.tx, 
                surf.tx = surf.tx),
  depth.bins = depth.bins,
  sim = x,
  sim.obs = obs
)

save(dive.sim, file = 'data/dive.sim.RData')
