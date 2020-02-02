data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

r.unif = max(outer(lambda, 2 * depth.bins$halfwidth, '/'))

m = dsdive.generator.matrix.uniformized(
  depth.bins = depth.bins, beta = beta, lambda = lambda, s0 = 2, 
  rate.uniformized = r.unif)

detach(dive.sim$params)
detach(dive.sim)