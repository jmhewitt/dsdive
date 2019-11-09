#' Gibbs sampler to estimate 3 stage dive model
#' 
#' Approximate the posterior distribution of parameters for the 3 stage dive 
#' model via Gibbs sampling over model parameters and imputed trajectories.
#' The dive bins and durations for the imputed trajectories are treated as 
#' fixed, but their stages are treated as latent, unknown variables that are 
#' sampled at each iteration if the stages are not provided.  
#' The full conditional posterior is used to update
#' the latent stages.  However, the model parameters are updated via adaptive
#' Metropolis-Hastings steps with a multivariate normal proposal distribution.
#' 
#' By default, the latent dive stages will be initialized by giving each dive 
#' stage equal time in the imputed trajectories.
#' 
#' @example examples/fit.fixedimputation.R
#' 
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
#' @param imputed.list A \code{list} of \code{dsdive} objects used as surrogate
#'   data for posterior approximation.
#' @param it number of iterations to run MCMC sampler for
#' @param burn number of iterations to discard at end of sampling
#' @param verbose If \code{verbose=TRUE}, then print diagnostic output while 
#'   sampling.
#' 
#' @export
#'
fit.fixedimputation = function(par, priors = rep(5, 12), imputed.list, it, 
                               burn = round(it/2), verbose = FALSE) {
  
  # update iterations to allow for initialization
  it = it+1
  
  # extract number of imputed dive trajectories
  N = length(imputed.list)
  
  # number of timepoints
  nt = length(imputed.list[[1]]$times)
  
  # initialize output
  res = list(
    stages = array(dim = c(it, N, nt)),
    beta = matrix(nrow = it, ncol = 6),
    lambda = matrix(nrow = it, ncol = 3),
    sub.tx = matrix(nrow = it, ncol = 2),
    surf.tx = matrix(nrow = it, ncol = 2)
  )
  
  # initialize parameters
  res$beta[1,] = par$beta
  res$lambda[1,] = par$lambda
  res$sub.tx[1,] = par$sub.tx
  res$surf.tx[1,] = par$surf.tx
  
  # initialize dive stages
  stage.vec = stagevec(length.out = nt, breaks = round(1:2 * nt/3))
  for(j in 1:N) {
    res$stages[1,j,] = stage.vec
  }
  
  # gibbs sample
  for(i in 2:it) {
    
    if(verbose) {
      message(paste('Starting iteration', i-1, sep = ' '))
    }
    
    # extract model parameters
    beta = matrix(res$beta[i-1,], nrow = 2)
    lambda = res$lambda[i-1,]
    sub.tx = res$sub.tx[i-1,]
    surf.tx = res$surf.tx[i-1,]
    
    # save model parameters
    res$beta[i,] = beta
    res$lambda[i,] = lambda
    res$sub.tx[i,] = sub.tx
    res$surf.tx[i,] = surf.tx
    
    # update stages for imputed trajectories
    for(j in 1:N) {
      
      if(verbose) {
        message(paste('  Updating stages for trajectory', j, sep = ' '))
      }
      
      # extract last stage vector
      stage.vec = res$stages[i-1, j, ]
      breaks = which(diff(stage.vec) == 1) + 1
      # update transition to second stage
      d = dsdive.ld.stages(breaks = breaks, fixed.ind =  2, beta = beta, 
                           lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
                           depths = imputed.list[[j]]$depths, 
                           durations = imputed.list[[j]]$durations, 
                           times = imputed.list[[j]]$times, 
                           depth.bins = depth.bins)
      breaks = c(sample(x = d$x, size = 1, prob = d$prob), breaks[2])
      stage.vec = stagevec(length.out = nt, breaks = breaks)
      # update transition to third stage
      d = dsdive.ld.stages(breaks = breaks, fixed.ind =  1, beta = beta, 
                           lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
                           depths = imputed.list[[j]]$depths, 
                           durations = imputed.list[[j]]$durations, 
                           times = imputed.list[[j]]$times, 
                           depth.bins = depth.bins)
      breaks = c(breaks[1], sample(x = d$x, size = 1, prob = d$prob))
      stage.vec = stagevec(length.out = nt, breaks = breaks)
      # save mcmc output
      res$stages[i, j, ] = stage.vec
    }
    
    if(verbose) {
      message('  Updating model parameters')
    }
    
    # Random walk update for model parameters
    
    # adapt random walk proposal distribution
    
    if(verbose) {
      message('  Adapting proposal distribution')
    }
  }
  
  # burn samples
  burn.inds = 1:burn
  res$stages = res$stages[-burn.inds,,, drop = FALSE]
  res$beta = res$beta[-burn.inds,]
  res$lambda = res$lambda[-burn.inds,]
  res$sub.tx = res$sub.tx[-burn.inds,]
  res$surf.tx = res$surf.tx[-burn.inds,]
  
  # return output
  res
}