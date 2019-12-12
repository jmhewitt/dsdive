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
#' @param priors List of parameters to specify prior distributions.  See 
#'   \code{dsdive.prior} for more details.
#' @param state.backup If not \code{NULL}, then a list that specifies a file 
#'   to which the sampler state will be dumped every \code{t} seconds.
#' @param scale.sigma.init Amount by which to scale the initial proposal 
#'   covariance matrices
#' 
#' @importFrom MHadaptive makePositiveDefinite
#' @importFrom stats optim rnorm runif var
#' 
#' @example examples/makeCompleteSingle.R
#' 
#' @export
#'
dsdive.fit.gibbs.cfg = function(cfg, it, verbose = FALSE, init, sigma,
                                priors, adapt = c(99, 100, 0.5, 0.75),
                                state.backup = list(t=Inf, file='state.RData'),
                                scale.sigma.init = 1) {
  
  # update number of iterations it to allow for initial parameters
  it = it + 1
  
  # build raw storage for posterior samples of model parameters
  par.init = params.toVec(par = init, spec = priors)
  n.par = c(2, 3) # 2 betas and three lambdas
  trace = matrix(nrow = it, ncol = sum(n.par))
  trace[1,] = par.init
  
  # store log-densities for data at each iteration
  ld = numeric(length = it)
  ld[1] = cfg$ld
  
  # compute cholesky decomposition for proposal covariance
  if(ncol(sigma[[1]]) + ncol(sigma[[2]]) != sum(n.par)) {
    stop(paste('The total parameter dimension', sum(n.par), 
               'differs from the proposal dimensions', 
               ncol(sigma[[1]]) + ncol(sigma[[2]]), 
               sep = ' '))
  } else {
    sigma.chol = list(
      t(chol(scale.sigma.init * sigma[[1]])),
      t(chol(scale.sigma.init * sigma[[2]]))
    )
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
  
  # store current log-prior value
  lp = dsdive.prior(par = init, spec = priors, log = TRUE)
  
  for(i in 2:it) {
    
    if(verbose) {
      tock = proc.time()
      message(paste('Iteration', i-1, sep=' '))
      message(paste('   step time:', round((tock - tick)[3],2), sep = ' '))
      tick = tock
    }
    
    # MH-RW sample depth bin transition parameters
    p = gibbs.mhrw.dsdive(x0 = trace[i-1,], ld0 = ld[i-1], lp0 = lp, 
                          sigma.chol = sigma.chol[[1]], priors = priors, 
                          ind = 1:2, verbose = verbose, cfg = cfg)
    trace[i,] = p$x
    ld[i] = p$ld
    lp = p$lp
    
    # MH-RW sample transition rates
    p = gibbs.mhrw.dsdive(x0 = trace[i,], ld0 = ld[i], lp0 = lp, 
                          sigma.chol = sigma.chol[[2]], priors = priors, 
                          ind = -(1:2), verbose = verbose, cfg = cfg)
    trace[i,] = p$x
    ld[i] = p$ld
    lp = p$lp
    
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
                                       var(trace[inds,1:2]))
      p.sigma2 = makePositiveDefinite((inds.len - 1) / inds.len * 
                                       var(trace[inds,-(1:2)]))
      # update proposal cholesky
      if(!(0 %in% p.sigma1) & !(0 %in% p.sigma2)) {
        sigma = list(p.sigma1, p.sigma2)
        sigma.chol = list(t(chol(p.sigma1)), t(chol(p.sigma2)))
      } 
    }
    
    # update trajectories, as necessary
    cfg = dsdive_cdtlimpute(cfg = cfg, params = p$params, i = i)
    ld[i] = cfg$ld
    
    # dump state
    if(!is.null(state.backup)) {
      tock = proc.time()[3]
      if(tock - tick.dump > state.backup$t) {
        
        if(verbose) {
          message('--- Saving state ---')
        }
        
        # initialize backup
        tmp = list(
          par = trace,
          ld = ld,
          sigma = sigma
        )
        
        # export or save imputed trajectories
        tmp = dsdive_outputImputed(cfg = cfg, output = tmp, 
                                   save.time = save.time, 
                                   file = state.backup$file)
        
        # save backup
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

  # export or save imputed trajectories
  res = dsdive_outputImputed(cfg = cfg, output = res, save.time = save.time, 
                             file = state.backup$file)

  res
}