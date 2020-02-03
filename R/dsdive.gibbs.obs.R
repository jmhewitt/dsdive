#' Likelihood for completely observed dive trajectories
#'
#' @param depths Depth bin indices in which the trajectory was observed
#' @param durations Amount of time spent in each depth bin
#' @param times Times at which each depth bin was entered
#' @param stages Record of dive stages at each depth bin
#' @param beta \eqn{2 x 3} matrix in which each column contains the diving 
#'  preference and directional persistence parameters for the DIVING, SUBMERGED, 
#'  and SURFACING dive stages.
#' @param lambda length 3 vector that specifies the transition rate, 
#'   respectively in the DIVING, SUBMERGED, and SURFACING stages.
#' @param sub.tx length 2 vector that specifies the first depth bin at which 
#'   transitions to the SUBMERGED stage can occur and the probability that such 
#'   a transition occurs at the next depth transition
#' @param surf.tx parameter that specifies the probability the trajectory will 
#'   transition to the SURFACING stage at the next depth transition
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param t0.dive Time at which dive started
#' @param d0.last If the depth bin that proceeded the first depth bin in 
#'   \code{depths}.  If the trajectory to be analyzed was started at the 
#'   surface, then set \code{c0.last=NULL}.
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#' 
#' @example examples/dsdive.gibbs.obs.R
#' 
#' 
dsdive.gibbs.obs = function(
  dsobs.list, t.stages.list, beta.init, lambda.init, verbose = FALSE, maxit, 
  checkpoint.fn, checkpoint.interval = 3600, pi1.prior, pi2.prior, 
  lambda1.prior, lambda2.prior, lambda3.prior, tstep, depth.bins) {

  
  lambda.priors.list = list(lambda1.prior, lambda2.prior, lambda3.prior)
  beta.priors.list = list(pi1.prior, pi2.prior)
  
  
  #
  # initialize sampler state and output
  #
  
  theta = list(beta = beta.init, lambda = lambda.init)
  
  trace = matrix(NA, nrow = maxit, ncol = 5)
  colnames(trace) = c('pi1', 'pi2', 'lambda1', 'lambda2', 'lambda3')
  
  tick.checkpoint = proc.time()[3]
  
  P.raw = lapply(1:3, function(s) {
    dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta.init, 
                        lambda = lambda.init, s0 = s, tstep = tstep, 
                        include.raw = TRUE)
  })
  
  
  
  #
  # gibbs sample
  #
  
  for(it in 1:maxit) {
    
    if(verbose) {
      message(paste('Iteration:', it, sep = ' '))
      tick = proc.time()[3]
    }
    
    
    #
    # sample model parameters from full conditional posterior densities
    #
    
    # update stage 1 parameters
    theta = dsdive.obs.sampleparams(
      dsobs.list = dsobs.list, t.stages.list = t.stages.list, P.raw = P.raw, 
      s0 = 1, depth.bins = depth.bins, beta = theta$beta, lambda = theta$lambda, 
      lambda.priors.list = lambda.priors.list, 
      beta.priors.list = beta.priors.list)
    
    # update stage 1 tx matrix 
    P.raw[[1]] = dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta.init, 
                                     lambda = lambda.init, s0 = 1, 
                                     tstep = tstep, include.raw = TRUE)
    # update stage 2 parameters
    theta = dsdive.obs.sampleparams(
      dsobs.list = dsobs.list, t.stages.list = t.stages.list, P.raw = P.raw, 
      s0 = 2, depth.bins = depth.bins, beta = theta$beta, lambda = theta$lambda, 
      lambda.priors.list = lambda.priors.list, 
      beta.priors.list = beta.priors.list)
    
    # update stage 2 tx matrix 
    P.raw[[2]] = dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta.init, 
                                     lambda = lambda.init, s0 = 2, 
                                     tstep = tstep, include.raw = TRUE)
    
    # update stage 3 parameters
    theta = dsdive.obs.sampleparams(
      dsobs.list = dsobs.list, t.stages.list = t.stages.list, P.raw = P.raw, 
      s0 = 3, depth.bins = depth.bins, beta = theta$beta, lambda = theta$lambda, 
      lambda.priors.list = lambda.priors.list, 
      beta.priors.list = beta.priors.list)
    
    # update stage 3 tx matrix 
    P.raw[[3]] = dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta.init, 
                                     lambda = lambda.init, s0 = 3, 
                                     tstep = tstep, include.raw = TRUE)
    
    if(verbose) {
      print(theta)
    }
    
    
    #
    # save trace
    #
    
    # save trace 
    trace[it,1:5] = unlist(theta[1:2])
    
    tock = proc.time()[3]
    
    if(verbose) {
      message(paste('   ', round(tock - tick, 1), 'seconds', sep = ' '))
    }
    
    
    #
    # checkpoint behaviors
    #
    
    # dump state
    if(tock - tick.checkpoint > checkpoint.interval) {
      if(verbose) {
        message('--- Checkpoint ---')
      }
      checkpoint.fn(trace = trace[1:it,], imputed.trace = imputed.trace[1:it])
      tick.checkpoint = proc.time()[3]
    }
    
  }
  
  list(
    theta = trace,
    dives = imputed.trace
  )
}