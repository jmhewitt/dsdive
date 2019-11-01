#' Builds a tidy data frame for observations of a dive
#' 
#' Converts raw dive information to a tidy dataframe that
#' is suitable for plotting with \code{ggplot2}.
#' 
#' @param depths Record of observed depth bins
#' @param times Times (in seconds) at which depth bin observations were made
#' @param stages Vector of guesses for which dive stage the trajectory was 
#'   in at each observation
#' @param depth.bins Vector that defines the depth bins (i.e., bin labels)
#' 
#' @importFrom stringr str_extract
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
  
  # convert seconds to minutes
  df$min = df$times/60
  
  # extract bin midpoints and ranges
  raw.depths = str_extract_all(depth.bins, '[0-9]+')
  depths.formatted = data.frame(do.call('rbind', 
                                        lapply(raw.depths, function(r) {
                                          r = as.numeric(r)
                                          c(r, mean(r))
                                        })))
  colnames(depths.formatted) = c('depth.min', 'depth.max', 'depth.mid')
  
  # enrich with depth bin information
  df = cbind(df, depths.formatted[df$depths+1,])
 
  df 
}