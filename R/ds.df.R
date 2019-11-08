#' Builds a tidy data frame for dive objects
#' 
#' Converts raw dive information to a tidy dataframe that
#' is suitable for plotting with \code{ggplot2}.  Most importantly, the function
#' converts depth bin indices to actual depths and depth ranges.
#' 
#' @param depths Record of observed depth bin indices
#' @param times Times (in seconds) at which depth bin observations were made
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param stages Vector of (guesses) for which dive stage the trajectory was 
#'   in at each observation
#' @param durations Vector specifying the amount of time spent in each depth bin
#' 
#' 
#' @example examples/ds.df.R
#' 
#' @export
#' 
ds.df = function(depths, times, depth.bins, stages = NULL, durations = NULL) {
  
  # initialize output
  df = data.frame(depths = depths, times = times)
  
  # add dive stage information
  if(!is.null(stages)) {
    df$stages = stages
  } else {
    df$stages = rep(NA, nrow(df))
  }
  df$stages = factor(df$stages)
  
  # add duration information
  if(!is.null(durations)) {
    if(length(depths) != length(durations)) {
      df$durations = c(durations, NA)
    } else {
      df$durations = durations
    }
  } else {
    df$durations = rep(NA, nrow(df))
  }
  df$durations.min = df$durations/60
  
  # convert seconds to minutes
  df$min = df$times/60
  
  # extract bin midpoints and ranges
  depths.formatted = cbind(depth.min = depth.bins[,1] - depth.bins[,2],
                           depth.max = depth.bins[,1] + depth.bins[,2],
                           depth.mid = depth.bins[,1])
  
  # enrich with depth bin information
  df = cbind(df, depths.formatted[df$depths,])
 
  df 
}