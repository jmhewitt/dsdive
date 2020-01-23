context("dsdive.augment.trajectory.R")

test_that("dive features are preserved when changing stage transition times", {
 
  data('dive.sim')
  attach(dive.sim)
  
  # recompute dive stages and durations when we change t.stages
  augmented = dsdive.augment.trajectory(depths = sim$depths, times = sim$times, 
                                        t.stages = 60*c(20,40))
  
  # augmented dive has "more" transitions because the new stage breaks do not 
  # coincide with depth transitions
  expect_gt(length(augmented$depths), length(sim$depths))
  
  # dive length remains unchanged 
  expect_equal(sum(sim$durations[is.finite(sim$durations)]),
               sum(augmented$durations[is.finite(augmented$durations)]))
  
  # the sequence of depth bin transitions also remains unchanged
  deltas = diff(sim$depths)
  deltas.augmented = diff(augmented$depths)
  expect_equal(deltas, deltas.augmented[deltas.augmented!=0])
  
  detach(dive.sim)
  
})
