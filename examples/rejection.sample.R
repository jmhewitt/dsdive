# number of samples to draw
n = 1e4

unif.circle = rejection.sample(
  n = n, 
  log.f = function(x) { log(sum(x^2) < 1) + log(pi) }, 
  log.g = function(x) 0, 
  rg = function(n) runif(2, min = -1, max = 1), 
  a = pi/4)

apply(unif.circle, 1, mean)

plot(t(unif.circle), xlim = c(-1,1), ylim = c(-1,1), asp = 1, 
     xlab = expression(X[1]), ylab = expression(X[2]))

