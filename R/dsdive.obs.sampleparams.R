#' Full conditional distribution for dive stages given data and parameters
#'
#'   
#' @param P.raw current collection of transition probability matrices
#' @param s0 the stage for which updated parameters should be sampled
#' @param beta the current depth bin transition model parameters
#' @param lambda the current depth bin transition rate model parameters
#'
#' @example examples/dsdive.obs.sampleparams.R
#' 
#' @importFrom bisque jac.log jac.logit
#' 
#' @export
#'
dsdive.obs.sampleparams = function(
  dsobs.list, t.stages.list, P.raw, s0, depth.bins, beta, lambda,
  lambda.priors.list, beta.priors.list) {

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
    
    if(any(theta.build$beta %in% c(0,1))) {
      -Inf
    } else {
      # update depth bin transition probability matrices
      P.raw[[s0]] = dsdive.obstx.matrix(
        depth.bins = depth.bins, beta = theta.build$beta, 
        lambda = theta.build$lambda, s0 = s0, tstep = tstep, 
        include.raw = TRUE)
      # likelihood
      dsdive.obsld(dsobs.list = dsobs.list, t.stages.list = t.stages.list, 
                   P.raw = P.raw, s0 = s0, sf = s0) + 
      # prior + jacobian
      theta.build$P + theta.build$J
    }
  }
  
  
  #
  # compute gaussian approximations to the posterior
  #
  
  if(s0==1) {
    x0 = c(log(lambda[1]), qlogis(beta[1]))
  } else if(s0==2) {
    x0 = log(lambda[2])
  } else if(s0==3) {
    x0 = c(log(lambda[3]), qlogis(beta[2]))
  }
  
  g = gaussapprox(logf = lp, init = x0)
  
  
  #
  # propose and accept/reject
  #
  
  # propose new model parameters
  x = g$rgaussapprox(n = 1)
  
  # metropolis ratio
  lR = lp(x) - lp(x0) + g$dgaussapprox(x = x0, log = TRUE) - 
    g$dgaussapprox(x = x, log = TRUE)
  
  # accept/reject
  if(log(runif(1)) <= lR) {
    res = x
  } else {
    res = x0
  }
  
  # back-transform, and return
  build.params(theta = x)[c('beta','lambda')]
}