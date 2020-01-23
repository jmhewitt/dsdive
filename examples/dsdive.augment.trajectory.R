data('dive.sim')
attach(dive.sim)

# recompute dive stages and durations when we change t.stages
augmented = dsdive.augment.trajectory(depths = sim$depths, times = sim$times, 
                                      t.stages = 60*c(20,40))

# dive before changing t.stages
pl1 = plot(x = sim, depth.bins = depth.bins)

# dive after changing t.stages
pl2 = plot(x = augmented, depth.bins = depth.bins)

detach(dive.sim)
 