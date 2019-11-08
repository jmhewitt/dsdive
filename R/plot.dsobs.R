#' Plotting functions for observations of dive trajectories
#' 
#' Builds a \code{ggplot2} graph.
#'
#' @param x \code{dsobs} object with record of observed depth bins
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param stages Optional vector of for which dive stage the trajectory was 
#'   in at each observation
#' @param errorbars If \code{TRUE}, then the minimum and maximum depth for each
#'   bin will be included in the plot.
#' @param underlay additional \code{ggplot2} layers to place before the main 
#'   plot layers
#' @param imputed.alpha transparency value used for plotting imputed 
#'   trajectories, if provided
#' @param imputed.list A \code{list} of \code{dsdive} objects to be plotted 
#'   underneath the main \code{dsdive} trajectory \code{x}.  The intent is that 
#'   \code{imputed.list} will contain imputed trajectories.
#' @param time.as.POSIXct if \code{TRUE}, will convert plotting times to POSIXct
#' @param ... (currently unused) additional plotting parameters
#' 
#' @example examples/plot.dsobs.R
#' 
#' @import ggplot2 ggthemes
#' 
#' @export
#'
plot.dsobs = function(x, depth.bins, stages = NULL, errorbars = FALSE, 
                      underlay = NULL, imputed.alpha = .3, imputed.list = NULL,
                      time.as.POSIXct = FALSE, ...) {
  
  # convert x to a plottable object
  df = ds.df(depths = x$depths, times = x$times, depth.bins = depth.bins, 
             stages = stages, time.as.POSIXct = time.as.POSIXct)

  #
  # build aes for depth bin observations
  #
  
  # base points
  m = aes(x = t.start, y = depth.mid)
  
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
  pl = ggplot(df)
  
  # underlay optional plot layers
  if(!is.null(underlay)) {
    pl = pl + underlay
  }
  
  # underlay imputed trajectories, if provided
  if(!is.null(imputed.list)) {
    
    nobs = length(x$times)
    
    m2 = aes(xmin = t.start, xmax = t.end, ymin = depth.min, 
             ymax = depth.max, fill = stages)
    
    # munge single imputed trajectory into a list format
    if(class(imputed.list) == 'dsdive') {
      imputed.list = list(imputed.list)
    }
    
    # add each imputed trajectory to plot
    for(i in 1:length(imputed.list)) {
      # build data frame for imputed trajectory
      df.imputed = ds.df(depths = imputed.list[[i]]$depths, 
                         times = imputed.list[[i]]$times, 
                         depth.bins = depth.bins, 
                         stages = imputed.list[[i]]$stages, 
                         durations = imputed.list[[i]]$durations,
                         time.as.POSIXct = time.as.POSIXct)
      # enrich duration information
      if(is.na(df.imputed$durations[nrow(df.imputed)])) {
          df.imputed$durations[nrow(df.imputed)] = 
            x$times[nobs] - df.imputed$times[nrow(df.imputed)]
          df.imputed$durations.min[nrow(df.imputed)] = 
            df.imputed$durations[nrow(df.imputed)] / 60
      }
      # add plot layer
      pl = pl + 
        geom_rect(mapping = m2, data = df.imputed, alpha = imputed.alpha)
    }
  }
  
  # main plot elements
  pl = pl + 
    # observations and depth bin ranges
    geom_pointrange(mapping = m, pch = 18) +
    # formatting
    scale_color_brewer('Dive stage', type = 'qual', palette = 'Set2') + 
    scale_fill_brewer('Dive stage', type = 'qual', palette = 'Set2') + 
    scale_y_reverse() + 
    xlab('Time') + 
    ylab('Depth (m)') + 
    theme_few() + 
    theme(panel.border = element_blank(),
          legend.position = 'top')
  
  pl
}