data('dive.sim')
attach(dive.sim)

# build partition of sim$times that 1) includes sim$times, and 2) has a 
# maximum increment of 100
times.refined = refine.partition(breaks = sim$times, max.width = 100)

# find where the original times appear in the refined partition
refined.mapping = sapply(sim$times, function(x) which(x==times.refined))

# display the difference in the partition
times.compare = rep(NA, length(times.refined))
times.compare[refined.mapping] = sim$times
cbind(original = times.compare, 
      refined = times.refined, 
      width = c(diff(times.refined), NA))

detach(dive.sim)
