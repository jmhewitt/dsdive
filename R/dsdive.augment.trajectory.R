#' Compute durations and stages from a complete record of depths and times
#' 
#' This function associates durations and stage labels with a complete 
#' observation record of all depth bin transitions, and stage transition times.
#' 
#' In practice, this is a support function that will facilitate re-computing 
#' durations given new stage transition times.
#' 
#' @example examples/dsdive.augment.trajectory.R
#' 
#' 
dsdive.augment.trajectory = function(depths, times, t.stages) {
  
  # there will be a depth bin and or stage transition at each of times.full
  times.full = sort(unique(c(times,t.stages)))

  #
  # package results
  #
  
  res = list(depths = depths[findInterval(times.full, times)], 
             times = times.full,
             durations = c(diff(times.full), NA), 
             stages = findInterval(times.full, t.stages) + 1)
  
  class(res) = 'dsdive'
  
  res
}