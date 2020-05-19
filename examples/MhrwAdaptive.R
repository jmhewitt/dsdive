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


# 
# Markov chain for two independent log-gamma r.v.'s
#

sampler = MhrwAdaptive$new(x = c(0, 0), 
                           mu = c(-1.27, -1.27), 
                           Sigma = diag(2),
                           lambda = c(1, 1), 
                           lp = function(x) sum(dlgamma(x = x, alpha = alpha, 
                                                    beta = beta, log = TRUE)))

x.samples = matrix(0, nrow = 1e4, ncol = 2)

for(i in 1:nrow(x.samples)) {
  x.samples[i,] = sampler$sample()$x
}

plot(density(x.samples[,1]))
curve(dlgamma(x = x, alpha = alpha, beta = beta), add = TRUE, col = 2)

sampler$print()
