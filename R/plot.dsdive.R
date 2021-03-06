#' Plotting functions for completely observed dive trajectories
#' 
#' Builds a \code{ggplot2} graph.
#'
#' @param x \code{dsdive} object with complete record of dive trajectory
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param imputed.alpha transparency value used for plotting imputed 
#'   trajectories, if provided
#' @param dsobs \code{dsobs} object with record of observed depth bins to 
#'   overlay on the main \code{dsdive} trajectory \code{x}
#' @param imputed.list A \code{list} of \code{dsdive} objects to be plotted 
#'   underneath the main \code{dsdive} trajectory \code{x}.  The intent is that 
#'   \code{imputed.list} will contain imputed trajectories.
#' @param underlay additional \code{ggplot2} layers to place before the main 
#'   plot layers
#' @param time.as.POSIXct if \code{TRUE}, will convert plotting times to POSIXct
#' @param ... (currently unused) additional plotting parameters
#' 
#' @example examples/plot.dsdive.R
#' 
#' @import ggplot2 ggthemes
#' 
#' @export
#'
plot.dsdive = function(x, depth.bins, imputed.alpha = .3, dsobs = NULL, 
                       imputed.list = NULL, underlay = NULL, 
                       time.as.POSIXct = FALSE, ...) {
  
  if(!inherits(x, 'dsdive')) {
    stop('x is not a dsdive object')
  }
  
  if(any(is.null(x$depths), is.null(x$times), is.null(x$stages))) {
    stop('x is not a valid dsdive object; missing depths, times, or stages')
  }

  # convert x to a plottable object
  df = ds.df(depths = x$depths, times = x$times, depth.bins = depth.bins, 
             stages = x$stages, durations = x$durations, 
             time.as.POSIXct = time.as.POSIXct)
  
  # enrich duration information
  if(is.na(df$durations[nrow(df)])) {
    if(!is.null(dsobs)) {
      nobs = length(dsobs$times)
      ntx = length(x$times)
      df$durations[nrow(df)] = dsobs$times[nobs] - x$times[ntx]
      df$durations.min[nrow(df)] = df$durations[nrow(df)] / 60
    }
  }
  
  
  #
  # build aes for complete depth bin records
  #
  
  # base elements
  m = aes(xmin = t.start, xmax = t.end, ymin = depth.min, 
          ymax = depth.max, fill = stages)
  
  # initialize plot
  pl = ggplot(df)
  
  # underlay optional plot layers
  if(!is.null(underlay)) {
    pl = pl + underlay
  }
  
  # underlay imputed trajectories, if provided
  if(!is.null(imputed.list)) {
    
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
        if(!is.null(dsobs)) {
          df.imputed$durations[nrow(df.imputed)] = 
            dsobs$times[nobs] - x$times[ntx]
          df.imputed$durations.min[nrow(df.imputed)] = 
            df.imputed$durations[nrow(df.imputed)] / 60
        }
      }
      # add plot layer
      pl = pl + geom_rect(mapping = m, data = df.imputed, alpha = imputed.alpha)
    }
  }
  
  # build main plot
  pl = pl + 
    # observations and depth bin ranges
    geom_rect(mapping = m) + 
    # formatting
    scale_fill_brewer('Dive stage', type = 'qual', palette = 'Set2') + 
    scale_y_reverse() + 
    xlab('Time') + 
    ylab('Depth (m)') + 
    theme_few() + 
    theme(panel.border = element_blank(),
          legend.position = 'top')
  
  # overlay observations, if provided
  if(!is.null(dsobs)) {
    
    # convert dsobs to a plottable object
    df.obs = ds.df(depths = dsobs$depths, times = dsobs$times, 
                   depth.bins = depth.bins, time.as.POSIXct = time.as.POSIXct)
    
    # base mapping
    m.obs = aes(x = t.start, y = depth.mid)
    
    # overlay observations
    pl = pl + 
      geom_point(mapping = m.obs, data = df.obs, 
                 inherit.aes = FALSE, pch = 18)
  }
  
  pl
}