#' Initialize a computing environment for use with gibbs sampling
#' 
#'
#' @param dives Vector of file paths OR a list of dive data.
#'   \itemize{
#'     \item If \code{dives} is a vector, then it must contain to file paths 
#'       to csv files that contain dive trajectory information.  The csv must 
#'       have labeled columns, and include \code{depths} and \code{times}
#'       columns.
#'       By default, the first time value will be used as 
#'       \code{t0.dive} during estimation.
#'     \item If \code{dives} is a list, then each entry must be a list that 
#'       contains \code{depths}, \code{times}, and \code{t0.dive}.
#'   }
#' @param depth.bins Vector of file paths OR a list of depth bins associated
#'   with dive data.  The \code{depth.bins} parameter is used in a similar 
#'   fashion as the \code{dives} parameter is.  The individual \code{depth.bins}
#'   objects or files must include columns or objects for depth bin 
#'   \code{center} and \code{halfwidth} values.  As a usage note, if multiple 
#'   \code{dives} need to use the same depth bin information, then the depth 
#'   bin file paths or objects should be duplicated.
#' @param cl A SNOW cluster object specifying nodes where computations will take
#'   place, and to which dive and depth bin information should be distributed.
#' 
#' @importFrom snow clusterApply
#' 
#' @export
#' 
#' @example examples/makeImputedLocal.R
#' 
makeImputedBatchDistributed = function(dives, depth.bins, cl, init, priors, it,
                                       inflation.factor.lambda) {
  
  # find the first dive id for each node's collection of dives
  batch.first = unique(ceiling(
    seq(from = 1, to = length(dives) + 1, length.out = length(cl) + 1)
  ))
  
  # generate local storage ids for dives
  ids = paste('dive', 1:length(dives), sep='')
  cluster.ids = paste('node', 1:length(cl), sep='')
  
  # merge data to load in batches on each node
  pkg = vector('list', length(cl))
  for(i in 1:length(cl)) {
    pkg[[i]] = list(cluster.id = cluster.ids[[i]])
  }
  for(i in 1:(length(batch.first)-1)) {
    inds = batch.first[i]:(batch.first[i+1]-1)
    pkg[[i]] = c(pkg[[i]], 
                 list(dive = dives[inds], depth.bins = depth.bins[inds], 
                    id = ids[inds], init = init, priors = priors, it = it, 
                    inflation.factor.lambda = inflation.factor.lambda)
                 )
  }
  
  # send data to nodes; set up local computing environments
  clusterApply(cl = cl, x = pkg, fun = function(p) {
    
    # ensure dsdive package is loaded
    require(dsdive)
    
    if(!is.null(p$dive)) {
      # initialize local computing environment for this dive
      cfg.local = makeImputedLocal(dives = p$dive, depth.bins = p$depth.bins, 
                                   init = p$init, priors = p$priors, it = p$it, 
                                   inflation.factor.lambda = 
                                     p$inflation.factor.lambda)
    } else {
      cfg.local = NULL
    }
      
    # save to worker's global environment
    assign(x = p$cluster.id, value = cfg.local, envir = globalenv())
    
  })
  
  # package configuration
  res = list(cluster.ids = cluster.ids, cl = cl)
  
  class(res) = 'dsImputedBatchDistributed'
  
  # compute initial ld
  params = params.toList(par = params.toVec(par = init, spec = priors), 
                         spec = priors)
  res$ld = dsdive_ld(cfg = res, params = params) + params$logJ
  
  res
}