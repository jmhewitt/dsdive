#' Split a sequence of integers into a partition with roughly equal-sized parts
#' 
#' Group the numbers \code{1:n.inds} into \code{n.partitions} such that all of 
#' the groups are about the same size.  This method can be useful for batching 
#' a complex operation into roughly equal-sized parts that will be evaluated in 
#' parallel.
#'   
#' @param n.inds Defines the number range \code{1:n.inds} that will be split
#'   into roughly equal-sized parts.
#' @param n.partitions The number of groups to create.
#'   
#' @example examples/balance_inds.R
#' 
#' @export
#' 
balance_inds = function(n.inds, n.partitions) {
  
  batch.size = floor(n.inds / n.partitions)
  remainder = n.inds - batch.size * n.partitions
  
  res = vector('list', n.partitions)
  
  start = 1
  for(i in 1:n.partitions) {
    end = start + batch.size + ifelse(i <= remainder, 1, 0) - 1
    if(end >= start) {
      res[[i]] = start:end
    }
    start = end + 1
  }
  
  res
}