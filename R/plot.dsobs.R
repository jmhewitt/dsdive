#' Plotting functions for observations of dive trajectories
#' 
#'
#' @param x \code{dsobs} object with record of observed depth bins
#' @param stages Optional vector of for which dive stage the trajectory was 
#'   in at each observation
#' @param depth.bins Vector that defines the depth bins (i.e., bin labels)
#' @param errorbars If \code{TRUE}, then the minimum and maximum depth for each
#'   bin will be included in the plot.
#' 
#' @example examples/plot.dsobs.R
#' 
#' @import ggplot2 ggthemes
#' 
#' @export
#'
plot.dsobs = function(x, depth.bins, stages = NULL, errorbars = FALSE, ...) {
  
  # convert x to a plottable object
  df = ds.df(depths = x$depths, times = x$times, depth.bins = depth.bins, 
             stages = stages)
  
  #
  # build aes for depth bin observations
  #
  
  # base points
  m = aes(x = min, y = depth.mid)
  
  # (optional) errorbars
  if(errorbars) {
    m$ymin = aes(depth.min)[[1]]
    m$ymax = aes(depth.max)[[1]]
  } else {
    m$ymin = aes(depth.mid)[[1]]
    m$ymax = aes(depth.mid)[[1]]
  }
  
  # colors for stages, if provided
  if(any(!is.na(df$stages))) {
    m$colour = aes(colour = stages)[[1]]
  }
  
  # build plot
  pl = ggplot(df) + 
    # observations and depth bin ranges
    geom_pointrange(mapping = m, pch = 18) +
    # formatting
    scale_color_brewer('Dive stage', type = 'qual', palette = 'Set2') + 
    scale_y_reverse() + 
    xlab('Time (min)') + 
    ylab('Depth (m)') + 
    theme_few() + 
    theme(panel.border = element_blank(),
          legend.position = 'top')
  
  pl
}