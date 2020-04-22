#' Initialize shared-memory pointers and node storage for Gibbs sampler
#' 
#' Initializes structures needed to run Gibbs sampler for covariate-based 
#' dive model across a shared-memory cluster.  Shared-memory clusters are 
#' important for allowing distributed computation of likelihoods with minimal 
#' communication costs and overhead.
#' 
#' @param cl Shared-memory cluster to be used to distribute some computations.
#'   The cluster must be fully initialized before running this function.
#'   The cluster requires random seeds and \code{Rdsm}-initialization.
#' @param envir Environment containing variables needed to send to remote nodes, 
#'   and used during initialization.  An error will be thrown if \code{envir} 
#'   does not contain all of the required variables.
#' 
#' @importFrom parallel clusterExport clusterEvalQ
#' @importFrom Rdsm mgrmakevar
#' 
#' @export
#' 
#' @example examples/dsdive.obs.sampleparams_shared.R
#'
gibbs_init_shared = function(cl, envir = environment()) {
  
  # data required by worker nodes
  cfg.static = c('dsobs.list', 'pi.designs', 'lambda.designs', 'beta1.prior', 
                 'beta2.prior', 'alpha1.prior', 'alpha2.prior', 'alpha3.prior', 
                 'tstep', 'depth.bins', 'T1.prior.params', 'T2.prior.params', 
                 'max.width', 'max.width.offset', 't0.prior.params', 
                 'tf.prior.params', 'delta')
  
  # data required for initialization
  cfg.init = c('t.stages.list', 'offsets', 'offsets.tf', 'beta.init', 
               'alpha.init')
  
  # verify data exists in envir
  missing = !(c(cfg.static, cfg.init) %in% names(envir))
  if(any(missing)) {
    stop(paste('envir does not contain required objects:', 
               paste(c(cfg.static, cfg.init)[missing], collapse = ', '),
               sep = ' '))
  }
    
  # verify list inputs are lists
  if(!inherits(envir$dsobs.list, 'list')) {
    stop('dsobs.list must be a list object')
  }
  if(!inherits(envir$t.stages.list, 'list')) {
    stop('dsobs.list must be a list object')
  }
  
  # send data to nodes
  clusterExport(cl, cfg.static, envir)
  
  
  #
  # allocate shared-memory objects; store access-copies in new environment
  #
  
  # new environment to store shared-memory objects
  shared = new.env()
  
  # node-level likelihoods
  do.call(what = mgrmakevar, envir = shared,
          args = list(varname = 'll.contrib', nr = length(cl), nc = 1,
                      vartype = 'double', savedesc = FALSE, cls = cl))
  
  # dive stage limits, for likelihood evaluation
  do.call(what = mgrmakevar, envir = shared,
          args = list(varname = 's0', nr = 1, nc = 1,
                      vartype = 'double', savedesc = FALSE, cls = cl))
  do.call(what = mgrmakevar, envir = shared,
          args = list(varname = 'sf', nr = 1, nc = 1,
                      vartype = 'double', savedesc = FALSE, cls = cl))
  
  # model parameters
  do.call(what = mgrmakevar, envir = shared,
          args = list(varname = 'beta1', nr = length(envir$beta.init[[1]]), 
                      nc = 1, vartype = 'double', savedesc = FALSE, cls = cl))
  do.call(what = mgrmakevar, envir = shared,
          args = list(varname = 'beta2', nr = length(envir$beta.init[[2]]), 
                      nc = 1, vartype = 'double', savedesc = FALSE, cls = cl))
  do.call(what = mgrmakevar, envir = shared,
          args = list(varname = 'alpha1', nr = length(envir$alpha.init[[1]]), 
                      nc = 1, vartype = 'double', savedesc = FALSE, cls = cl))
  do.call(what = mgrmakevar, envir = shared,
          args = list(varname = 'alpha2', nr = length(envir$alpha.init[[2]]), 
                      nc = 1, vartype = 'double', savedesc = FALSE, cls = cl))
  do.call(what = mgrmakevar, envir = shared,
          args = list(varname = 'alpha3', nr = length(envir$alpha.init[[3]]), 
                      nc = 1, vartype = 'double', savedesc = FALSE, cls = cl))
  
  # dive-specific random effects
  do.call(what = mgrmakevar, envir = shared,
          args = list(varname = 't.stages', nr = length(envir$t.stages.list), 
                      nc = 2, vartype = 'double', savedesc = FALSE, cls = cl))
  do.call(what = mgrmakevar, envir = shared,
          args = list(varname = 'offsets', nr = length(envir$t.stages.list), 
                      nc = 1, vartype = 'double', savedesc = FALSE, cls = cl))
  do.call(what = mgrmakevar, envir = shared,
          args = list(varname = 'offsets.tf', nr = length(envir$t.stages.list), 
                      nc = 1, vartype = 'double', savedesc = FALSE, cls = cl))
  
  
  #
  # initialize shared-memory values
  #
  
  # node-level likelihoods
  shared$ll.contrib[] = 0
  
  # dive stage limits, for likelihood evaluation
  shared$s0[] = 1
  shared$sf[] = 3
  
  # model parameters
  shared$beta1[] = envir$beta.init[[1]]
  shared$beta2[] = envir$beta.init[[2]]
  shared$alpha1[] = envir$alpha.init[[1]]
  shared$alpha2[] = envir$alpha.init[[2]]
  shared$alpha3[] = envir$alpha.init[[3]]
  
  # dive-specific random effects
  shared$t.stages[] = as.matrix(do.call(rbind, envir$t.stages.list))
  shared$offsets[] = envir$offsets
  shared$offsets.tf[] = envir$offsets.tf
  
  
  #
  # initialize state on nodes
  #
  
  clusterEvalQ(cl, {
    
    # enumerate dive id's this node is responsible for  
    local_ids = getidxs(length(dsobs.list))
    
    # number of dives for this node
    local_n = length(local_ids)
    
    # local storage for aligned dives
    dsobs.aligned = dsobs.list
    
    # initialize dive-specific transition matrices
    P.raw = lapply(1:local_n, function(local_ind) {
      lapply(1:3, function(s) {
        dsdive.obstxmat.cov(
          pi.designs = pi.designs, lambda.designs = lambda.designs, 
          beta1 = beta1[], beta2 = beta2[], alpha1 = alpha1[], 
          alpha2 = alpha2[], alpha3 = alpha3[], s0 = s, 
          ind = local_ids[local_ind], tstep = tstep, include.raw = TRUE, 
          depth.bins = depth.bins, delta = delta)
      })
    })
    
  })
  
  shared
}