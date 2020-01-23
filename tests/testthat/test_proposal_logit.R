context("Logit proposals")

test_that("Validate proposal distribution", {
  
  set.seed(2020)
  
  # define parameter support and initial conditions
  a = pi
  b = pi^3
  x0 = exp(2)
  
  # parameterize spread of proposal distribution
  sd = 1
  
  # number of monte carlo samples to draw
  m = 1e4
  
  # draw proposed x0 values
  x.samples = replicate(m, proposal.logit(x0 = x0, sd = sd, a = a, b = b)$x)
  
  # convert x0 to logit scale
  u = qlogis((x0-a)/(b-a))
    
  # # visually compare true proposal distribution to empirical samples
  # plot(density(qlogis((x.samples - a)/(b-a))))
  # curve(dnorm(x = x, mean = u, sd = sd), col = 2, add = TRUE,
  #       from = -5, to = 5)
  
  # verify proposals are centered around logit of x0
  rel_err = abs((mean(qlogis((x.samples-a)/(b-a))) - u)/u) * 100
  expect_lt(rel_err, 1)
})


test_that("Validate log-ratio", {
  
  set.seed(2020)
  
  # define parameter support and initial conditions
  a = pi
  b = pi^3
  x0 = exp(2)
  
  # parameterize spread of proposal distribution
  sd = 1
  
  # number of monte carlo samples to draw
  m = 1e5
  
  # draw proposed x0 values in Gibbs strategy; target dist'n. is uniform
  x.samples = numeric(m)
  x.samples[1] = 4
  for(i in 2:m) {
    # draw proposal
    prop = proposal.logit(x0 = x.samples[i-1], sd = sd, a = a, b = b)
    # accept/reject
    x.samples[i] = ifelse(log(runif(1)) < prop$lR, prop$x, x.samples[i-1])
  }
    
  # # visually look at density implied by samples
  # plot(density(x.samples))
  
  # verify mean and standard deviation
  mu = (a+b)/2
  v = (b-a)^2/12
  rel_err_mean = abs((mean(x.samples) - mu)/mu) * 100
  rel_err_var = abs((var(x.samples) - v)/v) * 100
  expect_lt(max(rel_err_mean, rel_err_var), 1)
})
