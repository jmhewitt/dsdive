#' Sample pi and lambda parameters for a single dive stage
#'
#' Sampler uses a Gaussian approximation to the full conditional posterior.
#'   
#' @param P.raw current collection of transition probability matrices
#' @param s0 the stage for which updated parameters should be sampled
#' @param beta the current depth bin transition model parameters
#' @param lambda the current depth bin transition rate model parameters
#' @param dsobs.list list of \code{dsobs} objects, which describe the 
#'   observation times and depths of a collection of dives
#' @param t.stages.list list of initial stage transition times for dives 
#'   observed in \code{dsobs.list}
#' @param P.raw list of continuous time probability transition matrices, and 
#'  components.
#' @param s0 stage for which model parameters should be sampled
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param beta directional preference model parameters.  
#'   See \code{dsdive.tx.params} for more details.
#' @param lambda diving rate model parameters.  See 
#'   \code{dsdive.tx.params} for more details.
#' @param lambda.priors.list list of parameters for the prior distributions for 
#'   the diving rate parameters.  See \code{dsdive.gibbs.obs} for more 
#'   information.
#' @param beta.priors.list list of parameters for the prior distributions for 
#'   the directional preference parameters.  See \code{dsdive.gibbs.obs} for  
#'   more information.  Note that these priors are specified as "pi" priors vs. 
#'   priors for "beta".
#' @param tstep Time between observations in \code{dsobs.unaligned}
#' @param gapprox (Optional) \code{gaussapprox} object containing the Gaussian 
#'   approximation used to propose model parameters.  If \code{NULL}, then a 
#'   Gaussian approximation will be computed.
#' @param output.gapprox \code{TRUE} to return the Gaussian approximation used 
#'   to propose model parameters.
#' @param delta If \code{delta>0}, then the probability transition matrices
#'   computed will use a transition matrix whose generator is 
#'   perturbed to allow much faster computation.  See \code{dsdive.obstx.matrix}
#'   for more details.
#'   
#' @example examples/dsdive.obs.sampleparams.R
#' 
#' @importFrom bisque jac.log jac.logit
#' @importFrom stats dgamma dbeta qlogis runif
#' 
#' @export
#'
dsdive.obs.sampleparams = function(
  dsobs.list, t.stages.list, P.raw, s0, depth.bins, beta, lambda,
  lambda.priors.list, beta.priors.list, tstep, gapprox = NULL,
  output.gapprox = FALSE, delta) {

  #
  # components for evaluating log-posteriors
  #
  
  build.params = function(theta) {
    
    #
    # assemble lambda parameters
    #
    
    # back-transform parameters
    L = exp(theta[1])
    # assemble parameters 
    lambda.eval = lambda
    lambda.eval[s0] = L
    # compute prior
    pl = lambda.priors.list[[s0]]
    P = dgamma(x = L, shape = pl[1], rate = pl[2], log = TRUE)
    # compute jacobian
    J = jac.log(x = L, log = TRUE)
    
    #
    # assemble beta parameters
    #
    
    if(is.na(theta[2])) {
      # assemble parameters
      beta.eval = beta
    } else {
      # back-transform parameters
      B = plogis(theta[2])
      # assemble parameters, and extract priors
      if(s0==1) {
        beta.eval = c(B, beta[2])
        pb = beta.priors.list[[1]]
      } else if(s0==3) {
        beta.eval = c(beta[1], B)
        pb = beta.priors.list[[2]]
      }
      # compute prior
      P = P + dbeta(x = B, shape1 = pb[1], shape2 = pb[2], log = TRUE)
      # compute jacobian
      J = J + jac.logit(x = B, log = TRUE)
    }
    
    # package results
    list(beta = beta.eval, lambda = lambda.eval, P = P, J = J)
  }
  
  
  lp = function(theta) {
    # Evaluate log-posterior at value of theta
    
    # extract parameters and related quantities
    theta.build = build.params(theta = theta)
     
    if(min(abs(outer(theta.build$beta, c(0,1), '-'))) <= 1e-16) {
      -Inf
    } else {
      # update depth bin transition probability matrices
      P.raw[[s0]] = dsdive.obstx.matrix(
        depth.bins = depth.bins, beta = theta.build$beta, 
        lambda = theta.build$lambda, s0 = s0, tstep = tstep, 
        include.raw = TRUE, delta = delta)
      # likelihood
      dsdive.obsld(dsobs.list = dsobs.list, t.stages.list = t.stages.list, 
                   P.raw = P.raw, s0 = s0, sf = s0) + 
      # prior + jacobian
      theta.build$P + theta.build$J
    }
  }
  
  
  #
  # compute or extract gaussian approximations to the posterior
  #
  
  if(s0==1) {
    x0 = c(log(lambda[1]), qlogis(beta[1]))
  } else if(s0==2) {
    x0 = log(lambda[2])
  } else if(s0==3) {
    x0 = c(log(lambda[3]), qlogis(beta[2]))
  }
  
  if(is.null(gapprox)) {
    g = gaussapprox(logf = lp, init = x0, method = 'Nelder-Mead')
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
  
  # back-transform parameters
  theta = build.params(theta = res)[c('beta','lambda')]
  
  
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
  
  res
}