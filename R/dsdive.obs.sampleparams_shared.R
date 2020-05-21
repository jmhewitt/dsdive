#' Sample pi and lambda parameters for a single dive stage
#'
#' Sampler uses a Gaussian approximation to the full conditional posterior.
#' Computations are carried out via shared-memory parallelization, and designed 
#' for use with the dive-specific covariate model.
#'   
#' @param s0 the stage for which updated parameters should be sampled
#' @param theta list containing current values of model parmaeters \code{beta1},
#'   \code{beta2}, \code{alpha1}, \code{alpha2}, and \code{alpha3}.
#' @param alpha.priors.list List of lists.  Each of the sublists specifies the 
#'   prior distribution for the parameter vectors \code{alpha1}, \code{alpha2}, 
#'   and \code{alpha3}.  See \code{dsdive.gibbs.obs.cov} for more details.
#' @param beta.priors.list List of lists.  Each of the sublists specifies the 
#'   prior distribution for the parameter vectors \code{beta1} and
#'   \code{beta2}.  See \code{dsdive.gibbs.obs.cov} for more details.
#' @param cl Shared-memory cluster to be used to distribute some computations.
#'   The cluster must be fully initialized before running this function.
#'   The cluster requires random seeds and \code{Rdsm}-initialization.
#' @param shared.env environment containing shared-memory variable pointers.  
#'   \code{shared.env} is expected to be the output from the 
#'   initialization function \code{gibbs_init_shared}.
#' @param gapprox (Optional) \code{gaussapprox} object containing the Gaussian 
#'   approximation used to propose model parameters.  If \code{NULL}, then a 
#'   Gaussian approximation will be computed.
#' @param output.gapprox \code{TRUE} to return the Gaussian approximation used 
#'   to propose model parameters.
#' @param optim.maxit maximum number of steps to take during numerical 
#'   optimization to compute Gaussian approximation to full conditional 
#'   posteriors used to propose model parameters
#' @param sample.betas \code{TRUE} to sample the directional preference 
#'   coefficients, or \code{FALSE} to sample the speed coefficients.
#' @param rw.sampler If \code{NULL}, then a new sampler will be created, 
#'   otherwise the sampler will be called.
#' @param adaptive \code{TRUE} to use adaptive Random walk Metropolis-Hastings
#'   proposals instead of Gaussian approximations
#' @param adaptation.frequency Random walk proposals will only be updated 
#'   at intervals of this step count
#' @param gapprox.approx_ld \code{TRUE} to use an approximate likelihood when 
#'   building Gaussian approximations to full conditional posteriors.
#'   
#' @example examples/dsdive.obs.sampleparams_shared.R
#' 
#' @importFrom stats dnorm runif
#' 
#' @export
#'
dsdive.obs.sampleparams_shared = function(
  s0, sample.betas, theta, alpha.priors.list, beta.priors.list, cl, 
  shared.env, gapprox = NULL, output.gapprox = FALSE, rw.sampler = NULL, 
  adaptive = FALSE, gapprox.approx_ld = FALSE, optim.maxit = 1e3, 
  adaptation.frequency = 10) {
  
  # build index subset for directional preference coefficients
  if(s0==1) {
    beta.inds = 1:length(theta$beta1)
  } else if(s0==2) {
    beta.inds = NULL
  } else if(s0==3) {
    beta.inds = 1:length(theta$beta2)
  }
  
  
  #
  # components for evaluating log-posteriors
  #
  
  build.params = function(theta_vec) {

    theta.eval = theta 
    
    # extract model parameters and compute prior
    if(s0==1) {
      if(sample.betas) {
        theta.eval$beta1 = as.numeric(theta_vec)
        P = sum(dnorm(x = theta.eval$beta1, mean = beta.priors.list[[1]]$mu, 
                      sd = beta.priors.list[[1]]$sd, log = TRUE))
      } else {
        theta.eval$alpha1 = as.numeric(theta_vec)
        P = sum(dnorm(x = theta.eval$alpha1, mean = alpha.priors.list[[1]]$mu, 
                      sd = alpha.priors.list[[1]]$sd, log = TRUE))
      }
    } else if(s0==2) {
      if(sample.betas) {
        stop('No beta coefficients to sample in stage 2!')
      }
      theta.eval$alpha2 = as.numeric(theta_vec)
      P = sum(dnorm(x = theta.eval$alpha2, mean = alpha.priors.list[[2]]$mu, 
                    sd = alpha.priors.list[[2]]$sd, log = TRUE))
    } else if(s0==3) {
      if(sample.betas) {
        theta.eval$beta2 = as.numeric(theta_vec)
        P = sum(dnorm(x = theta.eval$beta2, mean = beta.priors.list[[2]]$mu, 
                      sd = beta.priors.list[[2]]$sd, log = TRUE))
      } else {
        theta.eval$alpha3 = as.numeric(theta_vec)
        P = sum(dnorm(x = theta.eval$alpha3, mean = alpha.priors.list[[3]]$mu, 
                      sd = alpha.priors.list[[3]]$sd, log = TRUE))
      }
    }

    # package results
    list(theta = theta.eval, P = P)
  }
  
  
  lp = function(theta_vec) {
    # Evaluate log-posterior at value of theta_vec
    
    # extract parameters and related quantities
    theta.build = build.params(theta = theta_vec)
    
    # log-posterior
    dsdive.obsld_shared(theta = theta.build$theta, shared.env = shared.env, 
                        cl = cl, s0 = s0, sf = s0) + 
      theta.build$P
  }
  
  
  lp_approx = function(theta_vec) {
    # Evaluate approximate log-posterior at value of theta_vec
    
    # extract parameters and related quantities
    theta.build = build.params(theta = theta_vec)
    
    # log-posterior
    dsdive.obsld_approx_shared(theta = theta.build$theta, s0 = s0, sf = s0,
                               shared.env = shared.env, cl = cl) + 
      theta.build$P
  }
  
  
  
  #
  # compute or extract gaussian approximations to the posterior
  #
  
  if(s0==1) {
    if(sample.betas) {
      x0 = theta$beta1
    } else {
      x0 = theta$alpha1
    }
  } else if(s0==2) {
    x0 = theta$alpha2
  } else if(s0==3) {
    if(sample.betas) {
      x0 = theta$beta2
    } else {
      x0 = theta$alpha3
    }
  }
  
  if(adaptive) {
    
    if(is.null(rw.sampler)) {
      
      g = gaussapprox(logf = lp, init = x0, method = 'Nelder-Mead', 
                      control = list(fnscale =  -1, maxit = optim.maxit), 
                      optim.output = TRUE)
      
      rw.sampler = MhrwAdaptive$new(
        x = x0, 
        mu = g$optim.output$par, 
        Sigma = -solve(g$optim.output$hessian), 
        lambda = rep(1, length(x0)), lp = lp, C = .5, adaptive = TRUE, 
        adaptation_frequency = adaptation.frequency)
    }
    
    sampler_out = rw.sampler$sample()
    
    res = sampler_out$x
    accept = sampler_out$accepted
    
  } else {
    
    if(is.null(gapprox)) {
      if(gapprox.approx_ld) {
        g = gaussapprox(logf = lp_approx, init = x0, method = 'Nelder-Mead', 
                        control = list(fnscale =  -1, maxit = optim.maxit))
      } else {
        g = gaussapprox(logf = lp, init = x0, method = 'Nelder-Mead', 
                        control = list(fnscale =  -1, maxit = optim.maxit))
      }
    } else {
      g = gapprox
    }
    
    #
    # propose and accept/reject
    #
    
    # propose new model parameters
    x = g$rgaussapprox(n = 1)
    
    # metropolis ratio
    lR = lp(x) - lp(x0) + g$dgaussapprox(x = x0, log = TRUE) - 
      g$dgaussapprox(x = x, log = TRUE)
    
    # accept/reject
    accept = log(runif(1)) <= lR
    if(accept) {
      res = x
    } else {
      res = x0
    }
    
  }
  
  
  # back-transform parameters
  theta = build.params(theta_vec = res)$theta
  
  
  #
  # package results
  #
  
  res = list(
    theta = theta, 
    accepted = accept
  )
  
  if(output.gapprox) {
    res$g = g
  } 
  
  if(adaptive) {
    res$rw.sampler = rw.sampler
  }
  
  res
}