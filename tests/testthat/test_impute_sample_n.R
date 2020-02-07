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
  
})


test_that("sample num. uniformized depth bin tx's. across between-stage obs", {
  
  set.seed(2020)
  
  library(Matrix)
  library(dplyr)
  
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
  
  ind = max(which(sim.obs$times < t.stages[1]))
  
  d0 = sim.obs$depths[ind]
  df = sim.obs$depths[ind+1]
  
  s0 = sim.obs$stages[ind]
  sf = sim.obs$stages[ind+1]
  
  t0 = sim.obs$times[ind]
  tf = sim.obs$times[ind+1]
  
  # validate that selected transition is, indeed a within-stage transition
  expect_equal(s0+1, sf)
  
  # number of monte carlo samples to draw
  mcit = 1e4
  
  # MC sample number of uniformized depth bin tx's. between observations
  n.samples = replicate(n = mcit, expr = {
    n0 = dsdive.impute.sample_n(
      n0 = NULL, d0 = d0, df = df, s0 = s0, sf = sf, t0 = t0, tf = tf, 
      t.stages = t.stages, rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, 
      ff.s0 = NULL, ff.sf = NULL, n.bins = nrow(depth.bins), max.tx = 100, 
      ff.out = TRUE)
    
    n1 = dsdive.impute.sample_n(
      n0 = n0$n, d0 = d0, df = df, s0 = s0, sf = sf, t0 = t0, tf = tf, 
      t.stages = t.stages, rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, 
      ff.s0 = n0$ff.s0, ff.sf = n0$ff.sf, n.bins = nrow(depth.bins), 
      max.tx = 100)
    
    c(n0$n, n1)
  })
  
  # # MC sample number of uniformized depth bin tx's. between observations
  # n1.samples.cond = replicate(n = mcit, expr = {
  #   n0 = 4
  #   n1 = dsdive.impute.sample_n(
  #     n0 = n0, d0 = d0, df = df, s0 = s0, sf = sf, t0 = t0, tf = tf, 
  #     t.stages = t.stages, rate.unif = rate.unif, P.raw = P.raw, P.tx = P.tx, 
  #     ff.s0 = NULL, ff.sf = NULL, n.bins = nrow(depth.bins), 
  #     max.tx = 100)
  #   
  #   c(n0, n1)
  # })
  
  
  #
  # theoretical probability distributions, for validation
  #
  
  # depth bin distribution between t0 and stage transition time
  Pt.s0 = P.raw[[s0]]$evecs %*% 
          diag(exp(P.raw[[s0]]$evals * (t.stages[sf-1] - t0))) %*%
          P.raw[[s0]]$evecs.inv
  
  # depth bin distribution between stage transition time and tf
  Pt.sf = P.raw[[sf]]$evecs %*% 
    diag(exp(P.raw[[sf]]$evals * (tf - t.stages[sf-1]))) %*%
    P.raw[[sf]]$evecs.inv
  
  # depth bin distribution between t0 and tf
  Pt = Pt.s0 %*% Pt.sf
    
  # joint distribution
  f.joint = function(nmat) { apply(nmat, 1, function(n) {
    # Parameters: nmat is a matrix where each row holds c(n0, n1)
    
    # build n-step transition matrix for s0 transitions
    pmat.s0 = Diagonal(n = n.bins)
    if(n[1] > 0) {
      for(i in 1:n[1]) {
        pmat.s0 = pmat.s0 %*% P.tx[[s0]]
      }
    }
    
    # build n-step transition matrix for sf transitions
    pmat.sf = Diagonal(n = n.bins)
    if(n[2] > 0) {
      for(i in 1:n[2]) {
        pmat.sf = pmat.sf %*% P.tx[[sf]]
      }
    }
    
    # evaluate theoretical distribution at N = n
    (pmat.s0 %*% pmat.sf)[d0,df] * 
    dpois(x = n[2], lambda = rate.unif * (tf - t.stages[sf-1])) * 
    dpois(x = n[1], lambda = rate.unif * (t.stages[sf-1] - t0)) /
    Pt[d0,df]
  })}
  
  
  # MC estimates for joint pmf
  p.samples = table(n0 = n.samples[1,], n1 = n.samples[2,])/mcit
  n0.sampled = as.numeric(rownames(p.samples))
  n1.sampled = as.numeric(colnames(p.samples))
  
  # # MC estimates for conditional pmf
  # p.samples.cond = table(n0 = n1.samples.cond[1,], 
  #                        n1 = n1.samples.cond[2,])/mcit
  # n0.sampled.cond = as.numeric(rownames(p.samples.cond))
  # n1.sampled.cond = as.numeric(colnames(p.samples.cond))
  
  # theoretical joint distribution evaluated at observations
  joint.theoretical.support = expand.grid(n0 = n0.sampled, n1 = n1.sampled)
  joint.theoretical = f.joint(joint.theoretical.support)
  
  # helliger distance between MC estimates and theoretical joint dist'n.
  H.joint = sqrt(sum((sqrt(p.samples) - sqrt(joint.theoretical))^2)) / sqrt(2)
  
  # expect hellinger distance to be small
  expect_lt(H.joint, .1)
  
  # theoretical joint distribution evaluated at observations
  joint.theoretical.support.full = expand.grid(n0 = 0:40, n1 = 0:40)
  joint.theoretical.full = f.joint(joint.theoretical.support.full)

  # # visualize marginal distribution for n0
  # plot(n0.sampled, rowSums(p.samples),
  #      xlab = expression(n[0]), ylab = expression(P(N[0]==n[0])))
  # points(cbind(joint.theoretical.support.full,
  #              mass = joint.theoretical.full) %>%
  #          group_by(n0) %>% summarise(p = sum(mass)), col = 2)
  # 
  # # visualize marginal distribution for n1
  # plot(n1.sampled, colSums(p.samples),
  #      xlab = expression(n[1]), ylab = expression(P(N[1]==n[1])))
  # points(cbind(joint.theoretical.support.full,
  #              mass = joint.theoretical.full) %>%
  #          group_by(n1) %>% summarise(p = sum(mass)), col = 2)
  # 
  # # visualize all theoretical conditional distributions for n1 | n0
  # cbind(joint.theoretical.support.full, mass = joint.theoretical.full) %>%
  #   filter(n0 <=5) %>%
  #   mutate(n0 = factor(n0)) %>%
  #   group_by(n0) %>%
  #   mutate(p = mass/sum(mass)) %>%
  # ggplot(aes(x = n1, y = p, shape = n0, col = n0)) +
  #   geom_point()
  # 
  # # visualize a conditional distribution for n1 | n0
  # n0.cond = 0
  # pcond = p.samples[n0.sampled==n0.cond,]
  # pcond = pcond / sum(pcond)
  # plot(n1.sampled, pcond,
  #      xlab = expression(n[0]), ylab = expression(P(N[0]==n[0])))
  # points(cbind(joint.theoretical.support.full,
  #              mass = joint.theoretical.full) %>%
  #          filter(n0 == n0.cond) %>%
  #          mutate(p = mass/sum(mass)) %>% select(n1, p), col = 2)
  # 
  # # visualize an improved conditional distribution for n1 | n0
  # n0.cond = n0.sampled.cond
  # pcond = p.samples.cond[n0.sampled.cond==n0.cond,]
  # pcond = pcond / sum(pcond)
  # plot(n1.sampled.cond, pcond,
  #      xlab = expression(n[0]), ylab = expression(P(N[0]==n[0])))
  # points(cbind(joint.theoretical.support.full,
  #              mass = joint.theoretical.full) %>%
  #          filter(n0 == n0.cond) %>%
  #          mutate(p = mass/sum(mass)) %>% select(n1, p), col = 2)
  # 
  # # visualize a conditional distribution for n0 | n1
  # n1.cond = 10
  # pcond = p.samples[,n1.sampled==n1.cond]
  # pcond = pcond / sum(pcond)
  # plot(n0.sampled, pcond,
  #      xlab = expression(n[1]), ylab = expression(P(N[1]==n[1])))
  # points(cbind(joint.theoretical.support.full,
  #              mass = joint.theoretical.full) %>%
  #          filter(n1 == n1.cond) %>%
  #          mutate(p = mass/sum(mass)) %>% select(n0, p), col = 2)
  
  detach(params)
  detach(dive.sim)
  
})
