library(mvtnorm)

#
# univariate example
#

# specify distribution
mu = 1
sd = 3

# get gaussian approx. to... gaussian density
g = gaussapprox(logf = function(x) dnorm(x = x, mean = mu, sd = sd, log = TRUE), 
                init = 1, method = 'Brent', lower = -10, upper = 10, 
                optim.output = TRUE)

# sample from distribution
x = g$rgaussapprox(n = 1e4)

# see that mean and sd are recovered
c(mu = mean(x), sd = sd(x))

# density matches too
max(abs(g$dgaussapprox(x = x, log = FALSE) - 
          dnorm(x = x, mean = mu, sd = sd, log = FALSE)))


#
# multivariate example
#

# set dimension
p = 3

# generate random mean
mu = rnorm(p)

# generate random covariance matrix
sigma = cov(matrix(rnorm(5*p), ncol = p))

# get gaussian approx. to... gaussian density
g = gaussapprox(
  logf = function(x) dmvnorm(x = x, mean = mu, sigma = sigma, log = TRUE), 
  init = mu + 1, method = 'BFGS', optim.output = TRUE)

# sample from distribution
x = g$rgaussapprox(n = 1e5)

# see that mean and sd are recovered
max(abs(rowMeans(x) - mu))
max(abs((cov(t(x)) - sigma)))

# density matches too
max(abs(g$dgaussapprox(x = x, log = FALSE) - 
        dmvnorm(x = t(x), mean = mu, sigma = sigma, log = FALSE)))
