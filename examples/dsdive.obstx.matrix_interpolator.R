data('dive.sim')
attach(dive.sim)

# support for descent directional preferences
beta.seq = seq(.5, 1, length.out = 7)
beta.seq[length(beta.seq)] = .999

# support for descent speeds
lambda.seq = seq(.1, 2, length.out = 7)

# support for timesteps
tstep.seq = seq(0, 300, by = 100)

# build interpolating function
interpolator = dsdive.obstx.matrix_interpolator(
  depth.bins = depth.bins, beta.seq = beta.seq, lambda.seq = lambda.seq, 
  s0 = 1, tstep.seq = tstep.seq, m = 3, verbose = FALSE)

interpolator(beta = .8, lambda = 1, tstep = 300, i = 1, j = 4)

detach(dive.sim)
