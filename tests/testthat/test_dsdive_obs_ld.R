context("dsdive.obsld.R")

test_that("Verify implementation of exact likelihood function", {
  
  data('dive.sim')
  attach(dive.sim)
  attach(dive.sim$params)
  
  # extract time between observations
  tstep = diff(sim.obs$times[1:2])
  
  # get true stage transition times
  t.stages = sim$times[c(FALSE, diff(sim$stages)==1)]
  
  # build probability matrices for observations
  obstx.mat = lapply(1:3, function(s) {
    dsdive.obstx.matrix(depth.bins = depth.bins, 
                        beta = dive.sim$params$beta, 
                        lambda = dive.sim$params$lambda, 
                        s0 = s, tstep = tstep, include.raw = TRUE)
  })
  
  # compute likelihood of observations, given model parameters
  ld = dsdive.obsld(dsobs.list = list(sim.obs),
                    t.stages.list = list(t.stages), 
                    P.raw = obstx.mat, s0 = 1, sf = 3)
  
  expect_equal(ld, -17.6332728847878, tolerance = 1e-6)
  
  detach(dive.sim$params)
  detach(dive.sim)
  
})

test_that("Verify implementation of approximate likelihood function", {
  
  data('dive.sim')
  attach(dive.sim)
  attach(dive.sim$params)
  
  # extract time between observations
  tstep = diff(sim.obs$times[1:2])
  
  # get true stage transition times
  t.stages = sim$times[c(FALSE, diff(sim$stages)==1)]
  
  # build probability matrices for observations
  obstx.mat = lapply(1:3, function(s) {
    dsdive.obstx.matrix(depth.bins = depth.bins, 
                        beta = dive.sim$params$beta, 
                        lambda = dive.sim$params$lambda, 
                        s0 = s, tstep = tstep, include.raw = TRUE, 
                        delta = 1e-10)
  })
  
  # compute likelihood of observations, given model parameters
  ld = dsdive.obsld(dsobs.list = list(sim.obs),
                    t.stages.list = list(t.stages), 
                    P.raw = obstx.mat, s0 = 1, sf = 3)
  
  expect_equal(ld, -17.3078960893395, tolerance = 1e-6)
  
  detach(dive.sim$params)
  detach(dive.sim)
  
})
