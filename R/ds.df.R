#' Builds a tidy data frame for observations of a dive
#' 
#' Converts raw dive information to a tidy dataframe that
#' is suitable for plotting with \code{ggplot2}.
#' 
#' @param depths Record of observed depth bin indices
#' @param times Times (in seconds) at which depth bin observations were made
#' @param stages Vector of guesses for which dive stage the trajectory was 
#'   in at each observation
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' 
#' 
#' @example examples/ds.df.R
#' 
#' @export
#' 
ds.df = function(depths, times, depth.bins, stages = NULL) {
  
  # initialize output
  df = data.frame(depths = depths, times = times)
  
  # add dive stage information
  if(!is.null(stages)) {
    df$stages = stages
  } else {
    df$stages = rep(NA, nrow(df))
  }
  df$stages = factor(df$stages)
  
  # convert seconds to minutes
  df$min = df$times/60
  
  # extract bin midpoints and ranges
  depths.formatted = cbind(depth.min = depth.bins[,1] - depth.bins[,2],
                           depth.max = depth.bins[,1] + depth.bins[,2],
                           depth.mid = depth.bins[,1])
  
  # enrich with depth bin information
  df = cbind(df, depths.formatted[df$depths+1,])
 
  df 
}