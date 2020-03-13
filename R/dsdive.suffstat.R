#' Compute sufficient statistics for posterior distributions of complete dives
#'
#' @example examples/dsdive.suffstat.R
#' 
#' @export
#'   
dsdive.suffstat = function(depths, times, t.stages, depth.bins) {
  
  # get durations wrt. stages
  dive.augmented = dsdive.augment.trajectory(depths = depths, times = times, 
                                             t.stages = t.stages)
  
  # identify indices with valid durations
  finite.durations = is.finite(dive.augmented$durations)
  if(length(finite.durations)<length(dive.augmented$stages)) {
    finite.durations = c(finite.durations, FALSE)
  }
  
  # stage filters for transition rates
  s1 = dive.augmented$stages==1 & finite.durations
  s2 = dive.augmented$stages==2 & finite.durations
  s3 = dive.augmented$stages==3 & finite.durations
  
  # finite durations split by stages
  dur.stage.1 = dive.augmented$durations[s1]
  dur.stage.2 = dive.augmented$durations[s2]
  dur.stage.3 = dive.augmented$durations[s3]
  
  # depth bin widths split by stages
  width.stage.1 = 2 * depth.bins[dive.augmented$depths[s1], 2]
  width.stage.2 = 2 * depth.bins[dive.augmented$depths[s2], 2]
  width.stage.3 = 2 * depth.bins[dive.augmented$depths[s3], 2]
  
  
  # enumerate non-trivial, downward depth bin transitions by stage
  # output: (stage, downward?); all output is for non-trivial tx's.
  tx.nontrivial = sapply(2:length(dive.augmented$depths), function(i) {
    
    # get start and end depth bins
    d0 = dive.augmented$depths[i-1]
    df = dive.augmented$depths[i]
    
    # default output reflects no statistical information from transition
    res = c(NA,NA)
    
    # only update output if we observe a non-trivial depth bin transition
    if(d0!=df) {
      # get stage in which transition was made
      s = dive.augmented$stages[i]
      # get parameters for transition to determine if transition is trivial
      tx.params = dsdive.tx.params(depth.bins = depth.bins, d0 = d0, s0 = s, 
                                   beta = rep(.5,2), lambda = rep(1,3))
      # update output because transition is non-trivial
      if(length(tx.params$probs)>1) {
        res = c(s, df > d0)
      }
    }
    
    res
  })
  
  # identify depth bin transitions by stage
  s1.tx = tx.nontrivial[1,] == 1
  s3.tx = tx.nontrivial[1,] == 3
  
  # identify downward transitions by stage
  s1.tx.down = tx.nontrivial[2,s1.tx] == 1
  s3.tx.down = tx.nontrivial[2,s3.tx] == 1
  
  # summarize and package results
  list(
    # number of durations in each stage
    n.lambda = c(length(dur.stage.1), 
                 length(dur.stage.2), 
                 length(dur.stage.3)),
    # total normalized time spent in each stage
    d.lambda = c(sum(dur.stage.1 / width.stage.1),
                 sum(dur.stage.2 / width.stage.2),
                 sum(dur.stage.3 / width.stage.3)),
    # number of non-trivial transitions in each stage
    n.beta = c(sum(s1.tx, na.rm = TRUE), 
               sum(s3.tx, na.rm = TRUE)),
    # number of non-trivial down-transitions in each stage
    n.down = c(sum(s1.tx.down, na.rm = TRUE), 
               sum(s3.tx.down, na.rm = TRUE))
  )
}