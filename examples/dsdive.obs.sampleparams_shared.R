## 
## input for gibbs sampler
##

data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)
library(parallel)
library(Rdsm)

cl = makeCluster(2, 'SOCK')
clusterEvalQ(cl, library(dsdive))
clusterEvalQ(cl, library(Rdsm))

mgrinit(cl)

t.stages = sim$times[c(FALSE,diff(sim$stages)==1)]

tstep = diff(sim.obs$times[1:2])

obstx.mat = lapply(1:3, function(s) {
  dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta, 
                      lambda = lambda, s0 = s, tstep = tstep, 
                      include.raw = TRUE)
})

alpha.priors.list = list(
  list(mu = 0, sd = 3),
  list(mu = 0, sd = 3),
  list(mu = 0, sd = 3)
)

beta.priors.list = list(
  list(mu = 0, sd = 3),
  list(mu = 0, sd = 3)
)

T1.prior.params = c(25, .04)
T2.prior.params = c(56, .06)

dsobs.list = list(sim.obs, sim.obs)
t.stages.list = list(t.stages, t.stages)

beta.init = list(
  c(qlogis(.5), 0),
  c(0, 1)
)

alpha.init = c(beta.init, list(c(-1,-1)))

covs = data.frame(x1 = c(.5, 1), x2 = c(0, .3))

pi.formula = ~x1
lambda.formula = ~x1:x2


##
## barebones setup for gibbs sampler internals
##


n = length(dsobs.list)


#
# expand covariate design matrices
#

if(!inherits(pi.formula, 'list')) {
  pi.formula = list(pi.formula, pi.formula)
}

if(!inherits(lambda.formula, 'list')) {
  lambda.formula = list(lambda.formula, lambda.formula, lambda.formula)
}

pi.designs = lapply(pi.formula, function(f) model.matrix(f, covs))
lambda.designs = lapply(lambda.formula, function(f) model.matrix(f, covs))


#
# initialize nodes and shared memory pass-throughs
#

# required variables in environment
beta1.prior = beta.priors.list[[1]]
beta2.prior = beta.priors.list[[2]]
alpha1.prior = alpha.priors.list[[1]]
alpha2.prior = alpha.priors.list[[2]]
alpha3.prior = alpha.priors.list[[3]]
depth.bins = depth.bins
max.width = 100
max.width.offset = 30
t0.prior.params = c(1,1)
tf.prior.params = c(1,1)
delta = 1e-10
offsets = 0
offsets.tf = 0

shared.env = gibbs_init_shared(cl = cl, envir = environment())


#
# initialize sampler output
#

theta = list(beta1 = shared.env$beta1[], beta2 = shared.env$beta2[],
             alpha1 = shared.env$alpha1[], alpha2 = shared.env$alpha2[],
             alpha3 = shared.env$alpha3[])



##
## demo shared parameter sampling
##

# update stage 1 parameters
theta.raw = dsdive.obs.sampleparams_shared(
  s0 = 1, theta = theta, alpha.priors.list = alpha.priors.list,
  beta.priors.list = beta.priors.list, cl = cl, shared.env = shared.env)


##
## demo distributed likelihood evaluation
##

dsdive.obsld_shared(theta = theta, shared.env = shared.env, 
                    cl = cl, s0 = 1, sf = 3)


##
## clean up
##

stopCluster(cl)

detach(dive.sim$params)
detach(dive.sim)
