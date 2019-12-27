#' Compute sufficient statistics for posterior distributions of complete dives
#'
#'  @export
#'   
dsdive.suffstat = function(depths, durations, stages, depth.bins) {
  
  # identify indices with valid durations
  finite.durations = is.finite(durations)
  
  # stage filter for transitions
  s3.tx = stages==3
  
  # stage filters for transition rates
  s1 = stages==1
  s2 = stages==2 & finite.durations
  s3 = s3.tx & finite.durations
  
  # durations split by stages
  dur.stage.1 = durations[s1]
  dur.stage.2 = durations[s2]
  dur.stage.3 = durations[s3]
  
  # depths split by stages
  dep.stage.1 = depths[s1]
  dep.stage.2 = depths[s2]
  dep.stage.3 = depths[s3]
  dep.stage.3tx = depths[s3.tx]
  
  # depth bin widths split by stages
  width.stage.1 = 2 * depth.bins[dep.stage.1, 2]
  width.stage.2 = 2 * depth.bins[dep.stage.2, 2]
  width.stage.3 = 2 * depth.bins[dep.stage.3, 2]
  
  # max depth bin index
  n.bins = nrow(depth.bins)
  
  # identify non-trivial depth transitions
  nontrivial.tx1 = !(dep.stage.1 %in% c(1, n.bins))
  nontrivial.tx3 = !(dep.stage.3tx %in% c(1, n.bins))
  
  # tally non-trivial downward transitions
  n.down1 = sum(dep.stage.1[nontrivial.tx1] < 
                c(dep.stage.1, dep.stage.2[1])[c(FALSE, nontrivial.tx1)])
  tx.3 = dep.stage.3tx[nontrivial.tx3] < 
    c(dep.stage.3tx, NA)[c(FALSE, nontrivial.tx3)]
  n.down3 = sum(tx.3, na.rm = TRUE)
  
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
    n.beta = c(sum(nontrivial.tx1), sum(nontrivial.tx3) - sum(is.na(tx.3))),
    # number of non-trivial down-transitions in each stage
    n.down = c(n.down1, n.down3)
  )
}