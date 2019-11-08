#' Exploration of posterior surface via stochastic iteration
#'
#' Experimental method to explore the posterior distribution by approximating 
#' the posterior distribution by running MCMC over a family of imputed 
#' trajectories.  After the MCMC chain completes, the family of imputed 
#' trajectories is updated, and then the process is iterated until the poster
#' distribution approximation appears to converge wrt. Hellinger distance.
#' 
#' The method is experimental because its theoretical convergence properties 
#' are not fully explored.
#' 
#' @example examples/stochit.R
#' 
#' @importFrom MHadaptive Metro_Hastings
#' @importFrom bmk bmksensitive
#' 
#' @param depths Dive bins in which the trajectory was observed
#' @param t Times at which depths were observed
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param par List with parameter values at which to initialize the MCMC 
#'   algorithm. 
#'   \describe{
#'      \item{beta}{2x3 matrix with beta parameters}
#'      \item{lambda}{Vector of 3 elements}
#'      \item{sub.tx}{vector of depth index, and a probability}
#'      \item{surf.tx}{Vector of two parameters}
#'   }
#' @param priors Standard deviations for mean-0 normal priors on model 
#'   parameters, after transformation to real line
#' @param max.resample Maximum number of times to resample imputed trajectories
#' @param max.it Maximum number of MCMC iterations over parameters before 
#'   resampling imputed trajectories
#' @param tol Stop resampling imputed trajectories if Hellinger distance is 
#'   below \code{tol}
#' @param n.paths Number of trajectories to impute
#' @param n.par Number of parameter values to resample when imputing 
#'   trajectories
#' @param control.abc List with arguments of \code{dsdive.ldabc}, used for 
#'   imputing complete trajectories
#' @param control.mh List with arguments of \code{MHadaptive::Metro_Hastings}, 
#'   used for exploring the posterior distribution of model parameters, given 
#'   imputed trajectories
#' 
dsdive.stochit = function(depths, t, depth.bins, par, priors, max.resample, 
                          max.it, tol, n.paths, n.par, control.abc = NULL, 
                          control.mh = NULL, verbose = TRUE) {
  
  sub.tx1 = par$sub.tx[1]

  if(is.null(control.abc)) {
    control.abc = list(
      N = 1e1, 
      tries.max = 1e4,
      eps = 1,
      verbose = FALSE,
      steps.max = 1e3
    )
  }
  
  if(is.null(control.mh)) {
    control.mh = list(
      quiet = TRUE,
      burn_in = 1e3
    )
  }
  
  if(control.abc$N < n.paths + 1) {
    warning('Overriding default size of particle filter to support n.paths 
            argument.')
    control.abc$N = n.paths + 1
  }
  
  if(verbose) {
    message('Imputing initial family of trajectories')
  }
  
  # impute an initial family of trajectories
  paths.fam = dsdive.ldabc(beta = par$beta, lambda = par$lambda, 
                           sub.tx = par$sub.tx, surf.tx = par$surf.tx, 
                           depth.bins = depth.bins, 
                           steps.max = control.abc$steps.max, 
                           N = control.abc$N, depths = depths, t = t, 
                           tries.max = control.abc$tries.max, 
                           eps = control.abc$eps, dump.state = FALSE, 
                           verbose = control.abc$verbose, n.samples = n.paths, 
                           stage.init = 1)$sim
  
  # convenience function for evaluating log-posterior across the imputed 
  # family of trajectories
  lpost = function(theta) {
    # initialize log-posterior with prior information and offset
    res = sum(dnorm(x = theta, sd = priors, log = TRUE)) - n.paths*log(n.paths)
    
    # incorporate data
    for(i in 1:n.paths) {
      res = res + dsdive.ld(depths = paths.fam[[i]]$depths, 
                            durations = paths.fam[[i]]$durations, 
                            times = paths.fam[[i]]$times, 
                            stages = paths.fam[[i]]$stages, 
                            beta = matrix(theta[1:6], nrow = 2),
                            lambda = exp(theta[7:9]), 
                            sub.tx = c(sub.tx1, plogis(theta[10])), 
                            surf.tx = theta[11:12], 
                            depths.labels = depths.labels)
    }
    
    # return log-posterior
    res
  }
  
  if(verbose) {
    message('Computing initial approximation of posterior covariance matrix')
  }
  
  # optimize one of the paths in order to estimate hessian
  o = optim(par = c(par$beta, log(par$lambda), qlogis(par$sub.tx[2]), 
                    par$surf.tx), 
            fn = lpost, method = 'BFGS', hessian = TRUE, 
            control = list(fnscale = -1, reltol = .1))
  prop_sigma = -solve(o$hessian)
  
  # initialize structure for posterior comparison
  mh.last = list(trace = matrix(rep(o$par,2), nrow = 2))
  
  d = rep(tol, ncol(mh.last$trace))
  for(i in 1:max.resample) {
    
    if(verbose) {
      message(paste('Iteration:', i, sep = ' '))
      message(paste('  Initial Hellinger distances:', 
                    paste(round(d,2), collapse = ' '), sep = ' '))
    }
    
    # explore posterior distribution, conditional on imputed trajectories
    # Adaptive MH sampling of posterior for observed trajectory
    mh = Metro_Hastings(li_func = lpost, pars = o$par, prop_sigma = prop_sigma, 
                        iterations = max.it + control.mh$burn_in, 
                        quiet = control.mh$quiet, burn_in = control.mh$burn_in)
    
    # check convergence, resample trajectories if convergence not met
    d = bmksensitive(mh$trace, mh.last$trace)
    if(any(d>tol)) {
      
      if(verbose) {
        message(paste('  Updated Hellinger distances:', 
                      paste(round(d,2), collapse = ' '), sep = ' '))
        message('  Reimputing trajectories')
      }
      
      # select parameters to impute over
      inds = sample(x = nrow(mh$trace), size = n.paths, replace = FALSE)
      # update family of imputed trajectories
      for(j in 1:n.paths) {
        # extract parameters
        theta = mh$trace[inds[j],]
        # update an imputed trajectory
        paths.fam[[j]] = dsdive.ldabc(beta = matrix(theta[1:6], nrow = 2), 
                                      lambda = exp(theta[7:9]), 
                                      sub.tx = c(sub.tx1, plogis(theta[10])), 
                                      surf.tx = theta[11:12], 
                                      depths.labels = depths.labels, 
                                      steps.max = control.abc$steps.max, 
                                      N = control.abc$N, depths = depths, t = t, 
                                      tries.max = control.abc$tries.max, 
                                      eps = control.abc$eps, dump.state = FALSE, 
                                      verbose = control.abc$verbose, 
                                      n.samples = 1, stage.init = 1)$sim[[1]]
      }
      # swap posterior approximations
      mh.last = mh
    } else {
      break
    }
  }
  
  # return output from last posterior sampling
  mh
}