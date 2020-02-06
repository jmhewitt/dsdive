context("dsdive.impute.sample_n.R")

test_that("sample num. uniformized depth bin tx's. btwn within-stage obs", {
  
  set.seed(2020)
  
  library(Matrix)
  
  data('dive.sim')
  attach(dive.sim)
  attach(params)
  
  #
  # extract dive details
  #
  
  # number of depth bins
  n.bins = nrow(depth.bins)
  
  # time between observations
  tstep = diff(sim.obs$times[1:2])
  
  # stage transition times
  t.stages = sim$times[c(FALSE, diff(sim$stages)==1)]
  
  # uniformized transition rate
  rate.unif = max(outer(lambda, 2 * depth.bins[,2], '/'))
  
  # (continuous time) probability transition matrix for observations
  P.raw = lapply(1:3, function(s) {
    dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta, 
                        lambda = lambda, s0 = s, tstep = tstep, 
                        include.raw = TRUE)
  })
  
  # (discrete time) probability transition matrix for uniformized process
  P.tx = lapply(1:3, function(s) {
    dsdive.tx.matrix.uniformized(depth.bins = depth.bins, beta = beta, 
                                 lambda = lambda, s0 = s, 
                                 rate.uniformized = rate.unif)
  })
  
  
  #
  # within-stage transition example
  #
  
  ind = 2
  
  d0 = sim.obs$depths[ind]
  df = sim.obs$depths[ind+1]
  
  s0 = sim.obs$stages[ind]
  sf = sim.obs$stages[ind+1]
  
  t0 = sim.obs$times[ind]
  tf = sim.obs$times[ind+1]
  
  # validate that selected transition is, indeed a within-stage transition
  expect_equal(s0, sf)
  
  # number of monte carlo samples to draw
  mcit = 1e4
  
  # MC sample number of uniformized depth bin tx's. between observations
  n.samples = replicate(n = mcit, expr = {
    dsdive.impute.sample_n(
      n0 = NULL, d0 = d0, df = df, s0 = s0, sf = sf, t0 = t0, tf = tf, 
      t.stages = t.stages, rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, 
      ff.s0 = NULL, ff.sf = NULL, n.bins = nrow(depth.bins), max.tx = 100)
  })
  
  # theoretical probability distribution, for validation
  f = function(n) { sapply(n, function(n) {
    
    # build n-step transition matrix
    pmat = Diagonal(n = n.bins)
    for(i in 1:n) {
      pmat = pmat %*% P.tx[[s0]]
    }
    
    # evaluate theoretical distribution at N = n
    pmat[d0,df] * dpois(x = n, lambda = rate.unif * (tf - t0)) / 
      P.raw[[s0]]$obstx.mat[d0,df]
  })}
  

  # MC estimates for pmf
  p.samples = table(n.samples)/mcit
  
  # values of N for which MC estimates of pmf exist
  n.obs = as.numeric(names(p.samples))
  
  # evaluate theoretical distribution at values
  probs.theoretical = f(n.obs)
  
  # helliger distance between MC estimates and theoretical dist'n.
  H = sqrt(sum((sqrt(p.samples) - sqrt(probs.theoretical))^2)) / sqrt(2)
  
  # expect hellinger distance to be small
  expect_lt(H, .05)
  
  # # visualize differences between distributions
  # plot(n.obs, probs.theoretical)
  # points(n.obs, p.samples, col = 2, pch = 3)
  
  detach(params)
  detach(dive.sim)
  

  # 
  # #  
  # # use forward simulation to validate the target densities are correct
  # #
  # 
  # # N from forward-simulation of uniformized poisson process
  # n.samples.fwd = replicate(n = mcit, expr = {
  #   
  #   # rejection sample until we get a chain such that df.sample == df
  #   max.tries = 1e3
  #   for(it in 1:max.tries) {
  #     # draw number of transitions
  #     N = rpois(n = 1, lambda = rate.unif * (tf-t0))
  #     # initialize chain
  #     df.sample = d0
  #     # sample along chain
  #     if(N>0) {
  #       for(i in 1:N) {
  #         df.sample = sample(x = 1:n.bins, size = 1, 
  #                            prob = P.tx[[s0]][df.sample,])
  #       }
  #     }
  #     # stop sampling if we reach target df
  #     if(df.sample == df) {
  #       break
  #     }
  #   }
  #   
  #   ifelse(df.sample == df, N, NA)
  # })
  # 
  # table.samples.fwd = table(n.samples.fwd)
  # n.obs = as.numeric(names(table.samples.fwd))
  # probs.theoretical = f(n.obs)
  # probs.tgt = probs.theoretical > .05
  # max(abs(((table.samples.fwd/mcit - probs.theoretical)/probs.theoretical * 
  #            100)[probs.tgt]))
  # 
  # plot(n.obs, probs.theoretical)
  # points(n.obs, table.samples.fwd/mcit, col = 2, pch = 3)
  
  
  
})
