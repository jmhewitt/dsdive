#' Gibbs sampler for parameters of a model for dives across discrete depth bins
#'
#' @param dsobs.list list of \code{dsobs} objects, which describe the 
#'   observation times and depths of a collection of dives
#' @param t.stages.list list of initial stage transition times for dives 
#'   observed in \code{dsobs.list}
#' @param beta.init Initial values for directional preference model parameters.  
#'   See \code{dsdive.tx.params} for more details.
#' @param lambda.init Initial values for diving rate model parameters.  See 
#'   \code{dsdive.tx.params} for more details.
#' @param verbose \code{TRUE} to output sampler status while running
#' @param maxit number of Gibbs iterations to run
#' @param checkpoint.fn User-defined function to run during a checkpoint step;
#'   gives the user an opportunity to save partial output from the sampler
#' @param checkpoint.interval Number of seconds between calls to  
#'   \code{checkpoint.fn}
#' @param pi1.prior \code{shape1} and \code{shape2} parameters for Beta prior 
#'  distribution on the dive-stage model parameter \eqn{\pi^{(1)}}
#' @param pi2.prior \code{shape1} and \code{shape2} parameters for Beta prior 
#'  distribution on the ascent-stage model parameter \eqn{\pi^{(3)}}.  The 
#'  notation is a little odd because this is the SECOND preference parameter the 
#'  model estimates.
#' @param lambda1.prior \code{shape} and \code{rate} parameters for Gamma prior 
#'   distribution on the descent-stage diving rate \eqn{\lambda^{(1)}}.
#' @param lambda2.prior \code{shape} and \code{rate} parameters for Gamma prior 
#'   distribution on the bottom-stage diving rate \eqn{\lambda^{(2)}}.
#' @param lambda3.prior \code{shape} and \code{rate} parameters for Gamma prior 
#'   distribution on the ascent-stage diving rate \eqn{\lambda^{(3)}}.
#' @param tstep Time between observations in \code{dsobs.list}
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param T1.prior.params \code{shape} and \code{rate} parameters for Gamma 
#'   prior on the descent-stage duration.
#' @param T2.prior.params \code{shape} and \code{rate} parameters for Gamma 
#'   prior on the bottom-stage duration.
#' @param max.width The stage transition times are updated with a piecewise 
#'   proposal distribution.  \code{max.width} controls the maximum width of the 
#'   intervals for the proposal distribution.  This is a tuning parameter that 
#'   controls the numerical stability of the proposal distribution, which is 
#'   sampled via inverse CDF techniques.
#' @param max.width.offset The t0 and tf offsets are updated with a piecewise 
#'   proposal distribution.  \code{max.width.offset} controls the maximum width 
#'   of the intervals for the proposal distribution.  This is a tuning parameter 
#'   that controls the numerical stability of the proposal distribution, which 
#'   is sampled via inverse CDF techniques.
#' @param t0.prior.params \code{shape1} and \code{shape2} parameters for the 
#'   scaled and shifted Beta prior distribution for the t0 offset.
#' @param tf.prior.params \code{shape1} and \code{shape2} parameters for the 
#'   scaled and shifted Beta prior distribution for the tf offset.
#' @param offsets vector with initial values for t0 offsets.
#' @param offsets.tf vector with initial values for tf offsets.
#' @param warmup number of iterations during which the proposal distributions 
#'   will be updated at each step
#' @param cl cluster to be used to distribute some computations
#' @param delta If \code{delta>0}, then the probability transition matrices
#'   computed will use a transition matrix whose generator is 
#'   perturbed to allow much faster computation.  See \code{dsdive.obstx.matrix}
#'   for more details.
#'   
#' @example examples/dsdive.gibbs.obs.R
#' 
#' @importFrom parallel clusterApply
#' 
#' @export
#' 
dsdive.gibbs.obs = function(
  dsobs.list, t.stages.list, beta.init, lambda.init, verbose = FALSE, maxit, 
  checkpoint.fn, checkpoint.interval = 3600, pi1.prior, pi2.prior, 
  lambda1.prior, lambda2.prior, lambda3.prior, tstep, depth.bins, 
  T1.prior.params, T2.prior.params, max.width, max.width.offset,
  t0.prior.params, tf.prior.params, offsets, offsets.tf, warmup = Inf, 
  cl = NULL, delta) {
  
  n = length(dsobs.list)
  
  lambda.priors.list = list(lambda1.prior, lambda2.prior, lambda3.prior)
  beta.priors.list = list(pi1.prior, pi2.prior)
  
  # copy of dive observations with depths and timepoints aligned with t0=0
  dsobs.aligned = dsobs.list
  
  
  #
  # initialize sampler state and output
  #
  
  theta = list(beta = beta.init, lambda = lambda.init)

  trace = matrix(NA, nrow = maxit, ncol = 5)
  trace.t.stages = vector('list', length = maxit)
  trace.offsets = matrix(0, nrow = maxit, ncol = n)
  trace.offsets.tf = matrix(0, nrow = maxit, ncol = n)
  colnames(trace) = c('pi1', 'pi2', 'lambda1', 'lambda2', 'lambda3')
  
  tick.checkpoint = proc.time()[3]
  
  P.raw = lapply(1:3, function(s) {
    dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta.init, 
                        lambda = lambda.init, s0 = s, tstep = tstep, 
                        include.raw = TRUE, delta = delta)
  })
  
  # caches for proposal distributions
  proposaldists.theta = vector('list', 3)
  proposaldists.offsets = vector('list', n)
  proposaldists.offsets.tf = vector('list', n)
  
  
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
    theta.raw = dsdive.obs.sampleparams(
      dsobs.list = dsobs.aligned, t.stages.list = t.stages.list, P.raw = P.raw, 
      s0 = 1, depth.bins = depth.bins, beta = theta$beta, lambda = theta$lambda, 
      lambda.priors.list = lambda.priors.list, 
      beta.priors.list = beta.priors.list, tstep = tstep, 
      gapprox = proposaldists.theta[[1]], output.gapprox = TRUE, delta = delta)
    
    theta = theta.raw$theta

    # update stage 1 tx matrix 
    if(theta.raw$accepted) {
      P.raw[[1]] = dsdive.obstx.matrix(depth.bins = depth.bins, delta = delta,
                                       beta = theta$beta, 
                                       lambda = theta$lambda, s0 = 1, 
                                       tstep = tstep, include.raw = TRUE)
    }
    
    # # adapt stage 1 proposal distribution during warmup
    # if(it < warmup) {
    #   # only keep the proposal distribution if it generated an acceptance
    #   if(theta.raw$accepted) {
    #     proposaldists.theta[[1]] = theta.raw$g 
    #   } else {
    #     proposaldists.theta[1] = list(NULL)
    #   }
    # }
    # 
    # # lock in the final proposal distribution
    # if(it == warmup) {
    #   proposaldists.theta[[1]] = theta.raw$g
    # }
    
    # update stage 2 parameters
    theta.raw = dsdive.obs.sampleparams(
      dsobs.list = dsobs.aligned, t.stages.list = t.stages.list, P.raw = P.raw, 
      s0 = 2, depth.bins = depth.bins, beta = theta$beta, lambda = theta$lambda, 
      lambda.priors.list = lambda.priors.list, delta = delta,
      beta.priors.list = beta.priors.list, tstep = tstep,
      gapprox = proposaldists.theta[[2]], output.gapprox = TRUE)
    
    theta = theta.raw$theta
    
    # update stage 2 tx matrix 
    if(theta.raw$accepted) {
      P.raw[[2]] = dsdive.obstx.matrix(depth.bins = depth.bins, delta = delta,
                                       beta = theta$beta, 
                                       lambda = theta$lambda, s0 = 2, 
                                       tstep = tstep, include.raw = TRUE)
    }
    
    # # adapt stage 2 proposal distribution during warmup
    # if(it < warmup) {
    #   # only keep the proposal distribution if it generated an acceptance
    #   if(theta.raw$accepted) {
    #     proposaldists.theta[[2]] = theta.raw$g 
    #   } else {
    #     proposaldists.theta[2] = list(NULL)
    #   }
    # }
    # 
    # # lock in the final proposal distribution
    # if(it == warmup) {
    #   proposaldists.theta[[2]] = theta.raw$g
    # }
    
    
    # update stage 3 parameters
    theta.raw = dsdive.obs.sampleparams(
      dsobs.list = dsobs.aligned, t.stages.list = t.stages.list, P.raw = P.raw, 
      s0 = 3, depth.bins = depth.bins, beta = theta$beta, lambda = theta$lambda, 
      lambda.priors.list = lambda.priors.list, delta = delta,
      beta.priors.list = beta.priors.list, tstep = tstep,
      gapprox = proposaldists.theta[[3]], output.gapprox = TRUE)
    
    theta = theta.raw$theta
    
    # update stage 3 tx matrix 
    if(theta.raw$accepted) {
      P.raw[[3]] = dsdive.obstx.matrix(depth.bins = depth.bins, delta = delta,
                                       beta = theta$beta, 
                                       lambda = theta$lambda, s0 = 3, 
                                       tstep = tstep, include.raw = TRUE)
    }
    
    # # adapt stage 3 proposal distribution during warmup
    # if(it < warmup) {
    #   # only keep the proposal distribution if it generated an acceptance
    #   if(theta.raw$accepted) {
    #     proposaldists.theta[[3]] = theta.raw$g 
    #   } else {
    #     proposaldists.theta[3] = list(NULL)
    #   }
    # }
    # 
    # # lock in the final proposal distribution
    # if(it == warmup) {
    #   proposaldists.theta[[3]] = theta.raw$g
    # }
    
    if(verbose) {
      print(theta)
    }
    
    
    #
    # update stage transition time parameters and dive offsets
    #
    
    # # specify function to use to update dive-level random effects
    # if(is.null(cl)) {
    #   applyfn = lapply
    # } else {
    #   applyfn = function(X, FUN, ...) {
    #     clusterApply(cl = cl, x = X, fun = FUN, ...)
    #   }
    # }
    
    for(i in 1:n) {
      
      # sample dive stage transition times
      d = dsobs.aligned[[i]]
      t.stages.list[[i]] = dsdive.obs.sample.stages(
        depths = d$depths, times = d$times, t.stages = t.stages.list[[i]], 
        P.raw = P.raw, T.range = c(d$times[1], d$times[length(d$times)]), 
        depth.bins = depth.bins, T1.prior.params = T1.prior.params, 
        T2.prior.params = T2.prior.params, max.width = max.width, 
        debug = FALSE)$t.stages
      
      # sample dive start offset
      if(!is.null(t0.prior.params)) {
        
        d0 = dsdive.obs.sample.offsets(
          dsobs.aligned = d, dsobs.unaligned = dsobs.list[[i]], 
          offset = offsets[i], offset.tf = offsets.tf[i], 
          t.stages = t.stages.list[[i]], P.raw = P.raw, 
          depth.bins = depth.bins, tstep = tstep, max.width = max.width.offset, 
          prior.params = t0.prior.params, sample.start = TRUE, 
          lpapprox = proposaldists.offsets[[i]], output.lpapprox = TRUE)
        
        # extract new offset
        dsobs.aligned[[i]] = d0$dsobs.aligned
        offsets[i] = d0$offset
        
        # # adapt proposal distribution during warmup
        # if(it < warmup) {
        #   # only keep the proposal distribution if it generated an acceptance
        #   if(d0$accepted) {
        #     proposaldists.offsets[[i]] = d0$q1
        #   } else {
        #     proposaldists.offsets[i] = list(NULL)
        #   }
        # }
        # 
        # # lock in the final proposal distribution
        # if(it == warmup) {
        #   proposaldists.offsets[[i]] = d0$q1
        # }
        
      }
      
      # sample dive end offset
      if(!is.null(tf.prior.params)) {
        
        d0 = dsdive.obs.sample.offsets(
          dsobs.aligned = d, dsobs.unaligned = dsobs.list[[i]], 
          offset = offsets[i], offset.tf = offsets.tf[i], 
          t.stages = t.stages.list[[i]], P.raw = P.raw, 
          depth.bins = depth.bins, tstep = tstep, max.width = max.width.offset, 
          prior.params = tf.prior.params, sample.start = FALSE,
          lpapprox = proposaldists.offsets.tf[[i]], output.lpapprox = TRUE)
        
        # extract new offset
        dsobs.aligned[[i]] = d0$dsobs.aligned
        offsets.tf[i] = d0$offset
        
        # # adapt proposal distribution during warmup
        # if(it < warmup) {
        #   # only keep the proposal distribution if it generated an acceptance
        #   if(d0$accepted) {
        #     proposaldists.offsets.tf[[i]] = d0$q1
        #   } else {
        #     proposaldists.offsets.tf[i] = list(NULL)
        #   }
        # }
        # 
        # # lock in the final proposal distribution
        # if(it == warmup) {
        #   proposaldists.offsets.tf[[i]] = d0$q1
        # }
        
      }
      
    }
    
    
    #
    # save trace
    #
    
    trace[it,1:5] = unlist(theta[1:2])
    trace.t.stages[[it]] = t.stages.list
    trace.offsets[it,] = offsets
    trace.offsets.tf[it,] = offsets.tf
    
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
      checkpoint.fn(list(theta = trace[1:it,], 
                         trace.t.stages = trace.t.stages[1:it], 
                         trace.offsets = trace.offsets[1:it,],
                         trace.offsets.tf = trace.offsets.tf[1:it,]))
      tick.checkpoint = proc.time()[3]
    }
    
  }
  
  list(
    theta = trace,
    trace.t.stages = trace.t.stages,
    trace.offsets = trace.offsets,
    trace.offsets.tf = trace.offsets.tf
  )
}