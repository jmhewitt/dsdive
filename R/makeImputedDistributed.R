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
#' @example examples/makeImputedDistributed.R
#' 
makeImputedDistributed = function(dives, depth.bins, cl, init, priors, it,
                                   inflation.factor.lambda) {
  
  # generate local storage ids for dives
  ids = paste('dive', 1:length(dives), sep='')
  
  # merge data to send to nodes
  pkg = vector('list', length(dives))
  for(i in 1:length(dives)) {
    pkg[[i]] = list(dive = dives[[i]], depth.bins = depth.bins[[i]], 
                    id = ids[i], init = init, priors = priors, it = it, 
                    inflation.factor.lambda = inflation.factor.lambda)
  }
  
  # send data to nodes; set up local computing environments
  clusterApply(cl = cl, x = pkg, fun = function(p) {
    
    # ensure dsdive package is loaded
    require(dsdive)
    
    # if necessary, load dive and depth bin data from disk
    if(is.character(p$dive)) {
      # load dive and depth bin data from disk
      d = read.csv(file = p$dive, header = TRUE)
      db = read.csv(file = p$depth.bins, header = TRUE)
      # standardize column names for compatibility with dsdive conventions
      colnames(d) = tolower(colnames(d))
      colnames(db) = tolower(colnames(db))
      # convert data to list format
      p$dive = as.list(d)
      p$depth.bins = db
    }
    
    # initialize local computing environment for this dive
    cfg.local = makeImputedSingle(depth.bins = p$depth.bins, it = p$it, 
                                  depths = p$dive$depths, times = p$dive$times, 
                                  init = p$init, priors = p$priors, 
                                  inflation.factor.lambda = 
                                    p$inflation.factor.lambda, 
                                  t0.dive = p$dive$times[1], verbose = FALSE)
    
    # save to worker's global environment
    assign(x = p$id, value = cfg.local, envir = globalenv())
  })
  
  # package configuration
  res = list(cl = cl, ids = ids)
  
  class(res) = 'dsImputedDistributed'
  
  # compute initial ld
  params = params.toList(par = params.toVec(par = init, spec = priors), 
                         spec = priors)
  res$ld = dsdive_ld(cfg = res, params = params) + params$logJ
  
  res
}