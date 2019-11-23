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
#' @param sigma List of covariance matrices for block-Random walk proposals.
#'   First entry in list should be covariance matrix for depth transition 
#'   parameters, and second entry in list should be covariance matrix for 
#'   rate and stage transition parameters.
#' @param priors.sd vector of standard deviations for model parameters on their 
#'   respective transformed scales
#' @param state.backup If not \code{NULL}, then a list that specifies a file 
#'   to which the sampler state will be dumped every \code{t} seconds.
#' @param scale.sigma.init Amount by which to scale the initial proposal 
#'   covariance matrices
#' 
#' @importFrom MHadaptive makePositiveDefinite
#' 
#' @example examples/dsdive.fit.gibbs.R
#' 
#' @export
#'
dsdive.fit.gibbs = function(depths, times, durations = NULL, stages = NULL, 
                            depth.bins, t0.dive, it, verbose = FALSE, 
                            inflation.factor.lambda = 1.1, init, sigma = NULL,
                            priors.sd, 
                            adapt = c(100, 20, 0.5, 0.75),
                            state.backup = list(t=Inf, file='state.RData'),
                            scale.sigma.init = 1) {
  
  # update number of iterations it to allow for initial parameters
  it = it + 1
  
  # determine whether input represents a completely or partially observed dive
  partially.observed = is.null(durations)
  
  # build raw storage for posterior samples of model parameters
  par.init = params.toVec(init)
  n.par = c(3, 7)
  trace = matrix(nrow = it, ncol = sum(n.par))
  trace[1,] = par.init
  
  # store log-densities for data at each iteration
  ld = numeric(length = it)
  
  # compute initial jacobian
  logJ = params.toList(par = params.toVec(par = init))$logJ
  
  # if necessary, build storage for imputed trajectories, and initialize 
  #   latent trajectory
  if(partially.observed) {
    
    if(verbose) {
      message('Imputing initial trajectory')
    }
    
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
    trajectory$ld.true = trajectory$ld.true + logJ
    trajectory$ld = trajectory$ld + logJ
    
    # update log-densities
    ld[1] = trajectory$ld.true
    
    # save sample
    trace.imputed[[1]] = trajectory
    
  } else {
    # otherwise, set initial log density and trajectory object
    trajectory = list(
      depths = depths,
      stages = stages,
      times = times,
      durations = durations
    )
    
    # extract time at which stage 2 was entered, if any
    stage2.inds = which(trajectory$stages==2)
    if(length(stage2.inds)==0) {
      t.stage2tmp = NA
    } else {
      t.stage2tmp = trajectory$times[min(stage2.inds)]
    }
    
    ld[1] = dsdive.ld(depths = trajectory$depths, 
                      durations = trajectory$durations, 
                      times = trajectory$times, stages = trajectory$stages, 
                      beta = init$beta, lambda = init$lambda, 
                      sub.tx = init$sub.tx, surf.tx = init$surf.tx, 
                      depth.bins = depth.bins, t0.dive = t0.dive, 
                      t.stage2 = t.stage2tmp) + 
      logJ
  }
  
  
  # compute cholesky decomposition for proposal covariance
  if(!is.null(sigma)) { 
    # user-provided proposal covariance
    
    sigma.chol = list(
      t(chol(sigma[[1]])),
      t(chol(sigma[[2]]))
    )
    if(ncol(sigma.chol[[1]]) + ncol(sigma.chol[[2]]) != sum(n.par)) {
      stop(paste('The total parameter dimension', sum(n.par), 
                 'differs from the proposal dimensions', 
                 ncol(sigma.chol[[1]]) + ncol(sigma.chol[[2]]), 
                 sep = ' '))
    }
    
  } else { 
    # auto-tune proposal covariance based on trajectory if none given
    
    if(verbose) {
      message('Optimizing initial proposal covariance')
    }
    
    # extract time at which stage 2 was entered, if any
    stage2.inds = which(trajectory$stages==2)
    if(length(stage2.inds)==0) {
      t.stage2tmp = NA
    } else {
      t.stage2tmp = trajectory$times[min(stage2.inds)]
    }
    
    o = optim(par = trace[1,], fn = function(params) {
      # munge parameters
      params.prop.list = params.toList(par = params)
      # compute log-density
      prop.ld = dsdive.ld(depths = trajectory$depths,
                          durations = trajectory$durations,
                          times = trajectory$times, stages = trajectory$stages,
                          beta = params.prop.list$beta,
                          lambda = params.prop.list$lambda,
                          sub.tx = params.prop.list$sub.tx,
                          surf.tx = params.prop.list$surf.tx,
                          depth.bins = depth.bins, t0.dive = t0.dive, 
                          t.stage2 = t.stage2tmp) +
        params.prop.list$logJ
      # return log posterior
      prop.ld + sum(dnorm(x = params, sd = priors.sd, log = TRUE))
      
    }, method = 'BFGS', control = list(fnscale = -1), hessian = TRUE)
    
    sigma.tmp = -solve(o$hessian) * scale.sigma.init
    sigma = list(sigma.tmp[1:3,1:3], sigma.tmp[-(1:3),-(1:3)])
    sigma.chol = list(t(chol(sigma[[1]])), t(chol(sigma[[2]])))
    
  }
  
  if(verbose) {
    message('Sampling')
    tick = proc.time()
  }
  
  
  #
  # gibbs sample
  #
  
  if(!is.null(state.backup)) {
    tick.dump = proc.time()[3]
  }
  
  for(i in 2:it) {
    
    if(verbose) {
      tock = proc.time()
      message(paste('Iteration', i-1, sep=' '))
      message(paste('   step time:', round((tock - tick)[3],2), sep = ' '))
      tick = tock
    }
    
    # extract time at which stage 2 was entered, if any
    stage2.inds = which(trajectory$stages==2)
    if(length(stage2.inds)==0) {
      t.stage2tmp = NA
    } else {
      t.stage2tmp = trajectory$times[min(stage2.inds)]
    }
    
    # extract current model parameters
    params = params.toList(par = trace[i-1,])
    
    #
    # sample depth bin transition parameters
    #
    
    # propose new depth bin transition parameters
    params.prop = c(trace[i-1,1:3] + sigma.chol[[1]] %*% rnorm(n = 3),
                    trace[i-1,-(1:3)])
    params.prop.list = params.toList(par = params.prop)
    
    # compute proposed likelihood
    prop.ld = dsdive.ld(depths = trajectory$depths,
                        durations = trajectory$durations,
                        times = trajectory$times, stages = trajectory$stages,
                        beta = params.prop.list$beta,
                        lambda = params.prop.list$lambda,
                        sub.tx = params.prop.list$sub.tx,
                        surf.tx = params.prop.list$surf.tx,
                        depth.bins = depth.bins, t0.dive = t0.dive, 
                        t.stage2 = t.stage2tmp) +
      params.prop.list$logJ
    
    # accept/reject, and save parameters
    lR = prop.ld + sum(dnorm(x = params.prop, sd = priors.sd, log = TRUE)) - 
      (ld[i-1] + sum(dnorm(x = trace[i-1,], sd = priors.sd, log = TRUE)))
    
    if(log(runif(1)) <= lR) {
      if(verbose) {
        message(paste('   lR:', round(lR, 2), 'ACCEPTED', sep = ' '))
      }
      trace[i,] = params.prop
      ld[i] = prop.ld
      params = params.prop.list
    } else {
      if(verbose) {
        message(paste('   lR:', round(lR, 2), sep = ' '))
      }
      trace[i,] = trace[i-1,]
      ld[i] = ld[i-1]
    }
    
    
    #
    # sample transition rate and stage probability parameters
    #
    
    # propose new depth bin transition parameters
    params.prop = c(trace[i,1:3],
                    trace[i,-(1:3)] + sigma.chol[[2]] %*% rnorm(n = 7))
    params.prop.list = params.toList(par = params.prop)
    
    # compute proposed likelihood
    prop.ld = dsdive.ld(depths = trajectory$depths,
                        durations = trajectory$durations,
                        times = trajectory$times, stages = trajectory$stages,
                        beta = params.prop.list$beta,
                        lambda = params.prop.list$lambda,
                        sub.tx = params.prop.list$sub.tx,
                        surf.tx = params.prop.list$surf.tx,
                        depth.bins = depth.bins, t0.dive = t0.dive, 
                        t.stage2 = t.stage2tmp) +
      params.prop.list$logJ
    
    # accept/reject, and save parameters
    lR = prop.ld + sum(dnorm(x = params.prop, sd = priors.sd, log = TRUE)) - 
      (ld[i] + sum(dnorm(x = trace[i,], sd = priors.sd, log = TRUE)))
    
    if(log(runif(1)) <= lR) {
      if(verbose) {
        message(paste('   lR:', round(lR, 2), 'ACCEPTED', sep = ' '))
      }
      trace[i,] = params.prop
      ld[i] = prop.ld
      params = params.prop.list
    } else {
      if(verbose) {
        message(paste('   lR:', round(lR, 2), sep = ' '))
      }
      # no action, since leaving parameters untouched
    }
    
    
    # adapt RW proposal distributions
    if(i > adapt[1] && i %% adapt[2] == 0 && i < (adapt[4]*it) ) {   
      if(verbose) {
        message('   ADAPTING')
      }
      # select adaptation samples
      inds = floor(i*adapt[3]):i
      inds.len = length(inds)
      # get new proposal covariance
      p.sigma1 = makePositiveDefinite((inds.len - 1) / inds.len * 
                                       var(trace[inds,1:3]))
      p.sigma2 = makePositiveDefinite((inds.len - 1) / inds.len * 
                                       var(trace[inds,-(1:3)]))
      # update proposal cholesky
      if(!(0 %in% p.sigma1) & !(0 %in% p.sigma2)) {
        sigma = list(p.sigma1, p.sigma2)
        sigma.chol = list(t(chol(p.sigma1)), t(chol(p.sigma2)))
      } 
    }
    
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

      # get raw sampling weights
      W = exp(sapply(prop, function(p) p$w))
      # correct for edge cases in weights
      W[is.infinite(W)] = sign(W[is.infinite(W)])
      W[W==-1] = 0
      
      # sample dive 
      trajectory = prop[[sample(x = 2, size = 1, prob = W)]]
      
      # add jacobians
      trajectory$ld.true = trajectory$ld.true + params$logJ
      trajectory$ld = trajectory$ld + params$logJ
      
      # update log-density of proposal and trace
      trace.imputed[i] = list(trajectory)
      ld[i] = trajectory$ld.true
    }
    
    # dump state
    if(!is.null(state.backup)) {
      tock = proc.time()[3]
      if(tock - tick.dump > state.backup$t) {
        if(verbose) {
          message('--- Saving state ---')
        }
        tmp = list(
          par = trace,
          ld = ld,
          sigma = sigma
        )
        if(partially.observed) {
          tmp$trace.imputed = trace.imputed
        }
        save.time = date()
        save(tmp, save.time, file = state.backup$file)
        tick.dump = tock
      }
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
    ld = ld,
    sigma = sigma
  )
  
  if(partially.observed) {
    res$trace.imputed = trace.imputed[-1]
  }
  
  res
}