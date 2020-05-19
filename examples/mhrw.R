alpha = 2
beta = .5

# log-gamma denity
dlgamma = function(x, alpha, beta, log = FALSE) {
  res = (beta * x - exp(x) / alpha) - beta * log(alpha) - lgamma(beta)
  if(log) {
    res
  } else {
    exp(res)
  }
}

# single step MHRW
x = mhrw(x0 = 0, 
         lp = function(x) dlgamma(x = x, alpha = alpha, beta = beta, 
                                  log = TRUE), 
         cov.chol = 1)


# 
# Markov chain for two independent log-gamma r.v.'s
#

x.samples = matrix(0, nrow = 1e4, ncol = 2)

for(i in 2:nrow(x.samples)) {
  x.samples[i,] = mhrw(x0 = x.samples[i-1,], 
                       lp = function(x) sum(dlgamma(x = x, alpha = alpha, 
                                                    beta = beta, log = TRUE)),
                       cov.chol = diag(1,1))$x
}

plot(density(x.samples[,1]))
curve(dlgamma(x = x, alpha = alpha, beta = beta), add = TRUE, col = 2)


# 
# Markov chain for two independent cauchy r.v.'s
#

x.samples = matrix(0, nrow = 1e4, ncol = 2)

for(i in 2:nrow(x.samples)) {
  x.samples[i,] = mhrw(x0 = x.samples[i-1,], 
                       lp = function(x) sum(dcauchy(x = x, log = TRUE)),
                       cov.chol = diag(1,1))$x
}

plot(density(x.samples[,1]))
curve(dcauchy(x = x), add = TRUE, col = 2)
