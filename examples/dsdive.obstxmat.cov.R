data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

beta1 = c(0, 0)
beta2 = c(0, 0)
alpha1 = c(0, 0)
alpha2 = c(0, 0)
alpha3 = c(0, 0)

covs = data.frame(x1 = c(.5, 1), x2 = c(0, .3))

pi.formula = ~x1
lambda.formula = ~x1:x2


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
# build transition matrix
#

m = dsdive.obstxmat.cov(
  pi.designs = pi.designs, lambda.designs = lambda.designs, beta1 = beta1, 
  beta2 = beta2, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, s0 = 1, 
  ind = 1, tstep = 300, include.raw = TRUE, depth.bins = depth.bins, 
  delta = 1e-10)

detach(dive.sim$params)
detach(dive.sim)