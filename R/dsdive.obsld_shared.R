#' Likelihood for collection of partially observed dives
#' 
#' Computations are carried out via shared-memory parallelization, and designed 
#' for use with the dive-specific covariate model.
#' 
#' @param theta list containing current values of model parmaeters \code{beta1},
#'   \code{beta2}, \code{alpha1}, \code{alpha2}, and \code{alpha3}.
#' @param shared.env environment containing shared-memory variable pointers.  
#'   \code{shared.env} is expected to be the output from the 
#'   initialization function \code{gibbs_init_shared}.
#' @param cl Shared-memory cluster to be used to distribute some computations.
#'   The cluster must be fully initialized before running this function.
#'   The cluster requires random seeds and \code{Rdsm}-initialization.
#' @param s0 likelihood should include contributions for stages greater or 
#'   equal to \code{s0}
#' @param sf likelihood should include contributions for stages less than or 
#'   equal to \code{sf}
#'
#' @example examples/dsdive.obs.sampleparams_shared.R
#' 
#' @export
#'
dsdive.obsld_shared = function(theta, shared.env, cl, s0, sf) {
  
  # set stage limits for likelihood
  shared.env$s0[] = s0
  shared.env$sf[] = sf
  
  # push model parameters to nodes
  s.range = s0:sf
  if(1 %in% s.range) {
    shared.env$beta1[] = theta$beta1
    shared.env$alpha1[] = theta$alpha1
  }
  if(2  %in% s.range) {
    shared.env$alpha2[] = theta$alpha2
  }
  if(3 %in% s.range) {
    shared.env$beta2[] = theta$beta2
    shared.env$alpha3[] = theta$alpha3
  }
  
  clusterEvalQ(cl, {
    
    # local copy of likelihood range
    s.range = s0[]:sf[]
    
    # initialize local likelihood contribution
    ll = 0
    
    # loop over dives each node is responsible for
    for(local_ind in 1:local_n) {
      
      id = local_ids[local_ind]
      
      # update dive-specific transition matrices
      for(s in s.range) {
        P.raw[[local_ind]][[s]] = dsdive.obstxmat.cov(
          pi.designs = pi.designs, lambda.designs = lambda.designs, 
          beta1 = beta1[], beta2 = beta2[], alpha1 = alpha1[], 
          alpha2 = alpha2[], alpha3 = alpha3[], s0 = s, 
          ind = id, tstep = tstep, include.raw = TRUE, 
          depth.bins = depth.bins, delta = delta)
      }
      
      # dive-specific contribution to likelihood
      ll = ll + dsdive.obsld(dsobs.list = dsobs.list[id],
                             t.stages.list = list(t.stages[id,]),
                             P.raw = P.raw[[local_ind]], s0 = s0[], sf = sf[])
    }
    
    # push local likelihood contribution to main node
    ll.contrib[myinfo$id] = ll
  })
  
  # aggregate likelihood contributions from nodes
  sum(shared.env$ll.contrib[])
}
