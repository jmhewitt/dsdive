#' Merges repeated entries where we do not see a depth bin or stage transition
#' 
#' 
#' @example examples/dsdive.augment.trajectory.R
#' 
#' 
dsdive.simplify.trajectory = function(depths, times, stages) {
  
  # keep indices of places where we either change depth bins or stages
  inds.keep = c(TRUE, abs(diff(depths))==1) | c(TRUE, diff(stages)==1)
  
  #
  # package results
  #
  
  times.new = times[inds.keep]
  
  res = list(depths = depths[inds.keep], 
             times = times.new,
             durations = c(diff(times.new), NA),
             stages = stages[inds.keep])
  
  class(res) = 'dsdive'
  
  res
}