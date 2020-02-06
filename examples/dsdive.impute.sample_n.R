data('dive.sim')
attach(dive.sim)
attach(params)

#
# extract details about data
#

# time between observations
tstep = diff(sim.obs$times[1:2])

# stage transition indices
t.inds = which(c(FALSE, diff(sim$stages)==1))

# stage transition times
t.stages = sim$times[t.inds]

# uniformized transition rate
rate.unif = max(outer(lambda, 2 * depth.bins[,2], '/'))

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


#
# within-stage transition
#

ind = 2

d0 = sim.obs$depths[ind]
df = sim.obs$depths[ind+1]

s0 = sim.obs$stages[ind]
sf = sim.obs$stages[ind+1]

t0 = sim.obs$times[ind]
tf = sim.obs$times[ind+1]

document()
dsdive.impute.sample_n(
  n0 = NULL, d0 = d0, df = df, s0 = s0, sf = sf, t0 = t0, tf = tf, 
  t.stages = t.stages, rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, 
  ff.s0 = NULL, ff.sf = NULL, n.bins = nrow(depth.bins), max.tx = 100)


#
# between-stage transition
#

ind = max(which(sim.obs$times < t.stages[1]))

d0 = sim.obs$depths[ind]
df = sim.obs$depths[ind+1]

s0 = sim.obs$stages[ind]
sf = sim.obs$stages[ind+1]

t0 = sim.obs$times[ind]
tf = sim.obs$times[ind+1]

# sample number of depth bin transitions before stage 1->2 transition
document()
dsdive.impute.sample_n(
  n0 = NULL, d0 = d0, df = df, s0 = s0, sf = sf, t0 = t0, tf = tf, 
  t.stages = t.stages, rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, 
  ff.s0 = NULL, ff.sf = NULL, n.bins = nrow(depth.bins), max.tx = 100)


# sample number of depth bin transitions after stage 1->2 transition
document()
dsdive.impute.sample_n(
  n0 = 10, d0 = d0, df = df, s0 = s0, sf = sf, t0 = t0, tf = tf, 
  t.stages = t.stages, rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, 
  ff.s0 = NULL, ff.sf = NULL, n.bins = nrow(depth.bins), max.tx = 100)

detach(params)
detach(dive.sim)
