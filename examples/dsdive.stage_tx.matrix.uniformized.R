data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

r.unif = max(outer(lambda, 2 * depth.bins$halfwidth, '/'))

s0 = 1
sf = 3
s.range = s0:sf

# get standard uniformized transition matrices
tx.matrix = lapply(s.range, function(s) {
  dsdive.tx.matrix.uniformized(depth.bins = depth.bins, beta = beta, 
    lambda = lambda, s0 = s, rate.uniformized = r.unif)
})

m = dsdive.stage_tx.matrix.uniformized(m.list = tx.matrix)

detach(dive.sim$params)
detach(dive.sim)