#' Gibbs sampler to approximate posterior for model parameters
#' 
#' The code is designed to approximate posterior distributions when the dive 
#' trajectory is either completely or incompletely observed.  When the 
#' trajectory is incompletely observed, a latent dive trajectory will be updated 
#' using an independence proposal at each iteration.
#'
#' 
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param depths record of depth bins the trajectory should visit
#' @param times times at which the depth bins should be visited
#' @param s0 dive stage at which the trajectory should be started from
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
#' @param inflation.factor.lambda In order to facilitate bridged transitions, 
#'   the transition rate of the overall process must be inflated to allow the 
#'   possibility of self-transitions.  Self-transitions allow bridged paths to 
#'   dynamically modify the total number of transitions between observed values
#'   so that a valid path between observations is always possible.  The 
#'   \code{inflation.factor.lambda} parameter implicitly controls the number of 
#'   self-transitions that will occur.  Larger values will create more 
#'   self-transitions.
#' @param verbose If \code{TRUE}, then the sampler's progress will be printed 
#'   during sampling.
#' @param precompute.bridges If \code{TRUE}, then the bridged transition 
#'   matrices will be precomputed.  Enabling this option will increase the 
#'   memory overhead of the method, but will reduce its runtime.
#' @param t0.dive Time at which dive started
#' @param init List of parameter values at which to initialize the MCMC chain
#' @param sigma Covariance matrix for Random walk proposals
#' 
#' @example examples/dsdive.impute.R
#' 
#' @export
#'
dsdive.fit.gibbs = function(depths, times, durations = NULL, stages = NULL, 
                            depth.bins, t0.dive, it, verbose = FALSE, 
                            inflation.factor.lambda = 1.1, init, sigma) {
  
  # update number of iterations it to allow for initial parameters
  it = it + 1
  
  # determine whether input represents a completely or partially observed dive
  partially.observed = is.null(durations)
  
  # build raw storage for posterior samples of model parameters
  par.init = params.toVec(init)
  n.par = length(par.init)
  trace = matrix(nrow = it, ncol = n.par)
  trace[1,] = par.init
  
  # compute cholesky decomposition for proposal covariance
  sigma.chol = t(chol(sigma))
  if(ncol(sigma.chol) != n.par) {
    stop(paste('The parameter dimension', n.par, 
               'differs from the proposal dimension', ncol(sigma.chol), 
               sep = ' '))
  }
  
  # store log-densities for data at each iteration
  ld = numeric(length = it)
  
  # if necessary, build storage for imputed trajectories, and initialize 
  #   latent trajectory
  if(partially.observed) {
    
    # initialize record of latent trajectories
    trace.imputed = vector('list', it)
    
    # sample first latent trajectory
    trajectory = dsdive.fastimpute(M = 1, depth.bins = depth.bins, 
                                   depths = depths, times = times, 
                                   s0 = 1, beta = init$beta, 
                                   lambda = init$lambda, 
                                   sub.tx = init$sub.tx, 
                                   surf.tx = init$surf.tx, 
                                   inflation.factor.lambda = 
                                     inflation.factor.lambda, 
                                   verbose = FALSE, 
                                   precompute.bridges = TRUE, 
                                   t0.dive = t0.dive, 
                                   resample = FALSE)[[1]]
    
    # add log jacobians to sample
    logJ = params.toList(par = params.toVec(par = init), 
                         sub.tx1 = init$sub.tx[1])$logJ
    trajectory$ld.true = trajectory$ld.true + logJ
    trajectory$ld = trajectory$ld + logJ
    
    # update log-densities
    ld[1] = trajectory$ld.true
    
    # save sample
    trace.imputed[[1]] = trajectory
    
  } else {
    # otherwise, set initial log density and trajectory object
    
  }
  
    
  # gibbs sample
  for(i in 2:it) {
    
    if(verbose) {
      message(paste('Iteration', i-1, sep=' '))
    }
    
    # extract current model parameters
    params = params.toList(par = trace[i-1,], sub.tx1 = 1)
    
    # propose new model parameters
    # params.prop = trace[i-1,] + sigma.chol %*% rnorm(n = n.par)
    # params.prop.list = params.toList(par = params.prop, sub.tx1 = 1)
    # 
    # # compute proposed likelihood
    # prop.ld = dsdive.ld(depths = trajectory$depths, 
    #                     durations = trajectory$durations, 
    #                     times = trajectory$times, stages = trajectory$stages, 
    #                     beta = params.prop.list$beta, 
    #                     lambda = params.prop.list$lambda, 
    #                     sub.tx = params.prop.list$sub.tx, 
    #                     surf.tx = params.prop.list$surf.tx, 
    #                     depth.bins = depth.bins, t0.dive = t0.dive) + 
    #   params.prop.list$logJ
    
    # accept/reject
    
    
    # params.prop = params
    
    # don't forget to add jacobians!!!!
    
    # accept/reject
    # accept/reject
    # lR = prop$ld.true - prop$ld - (lds[i-1,1] - lds[i-1,2])
    # if(log(runif(1)) <= lR) {
    #   trace.imputed[[i]] = prop
    #   trajectory = prop
    # } else {
    #   trace.imputed[[i]] = trace.imputed[[i-1]]
    # }
    
    # save parameters
    trace[i,] = params.toVec(par = params)
    
    # adapt RW proposal distribution
    
    
    if(partially.observed) {
      # propose trajectory
      prop = dsdive.fastimpute(M = 2, depth.bins = depth.bins, depths = depths, 
                               times = times, s0 = 1, beta = params$beta, 
                               lambda = params$lambda, sub.tx = params$sub.tx, 
                               surf.tx = params$surf.tx, 
                               inflation.factor.lambda = 
                                 inflation.factor.lambda, verbose = FALSE, 
                               precompute.bridges = TRUE, 
                               t0.dive = t0.dive, 
                               trajectory.conditional = trajectory)

      # sample dive
      W = exp(sapply(prop, function(p) p$w))
      trajectory = prop[sample(x = 2, size = 1, prob = W)]
      
      # add jacobians
      trajectory$ld.true = trajectory$ld.true + params$logJ
      trajectory$ld = trajectory$ld + params$logJ
      
      # update log-density of proposal and trace
      trace.imputed[i] = trajectory
      ld[i] = trajectory$ld.true
    }
    
  }
  
  # remove initial parameter values
  trace = trace[-1,]
  ld = ld[-1]
  
  
  #
  # package results
  #
  
  res = list(
    par = trace,
    ld = ld
  )
  
  if(partially.observed) {
    res$trace.imputed = trace.imputed[-1]
  }
  
  res
}