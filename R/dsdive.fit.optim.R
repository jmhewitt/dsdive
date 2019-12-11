#' Find posterior mode and hessian for fixed configuration
#' 
#' Useful for revising initial parameters and proposal distributions for Gibbs 
#' sampler.
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
#' @importFrom stats optim
#' 
#' @example examples/dsdive.fit.optim.R
#' 
#' @export
#'
dsdive.fit.optim = function(cfg, method = 'BFGS', maxit = 1e3, verbose = FALSE, 
                            init, priors, reltol = sqrt(.Machine$double.eps),
                            hessian  = TRUE) {

  par.init = params.toVec(par = init, spec = priors)
  
  o = optim(par = par.init, fn = function(params) {
    # munge parameters
    params.prop.list = params.toList(par = params, spec = priors)
    # compute log-density
    prop.ld = dsdive_ld(cfg = cfg, params = params.prop.list) +
      params.prop.list$logJ
    # compute log posterior
    r = prop.ld + 
      dsdive.prior(par = params.prop.list, spec = priors, log = TRUE)
    # diagnostics
    if(verbose) {
      print(params)
      message(r)
    }
    # return
    r
  }, method = method, hessian = hessian,
  control = list(maxit = maxit, fnscale = -1, reltol = reltol))
  
  r = list(par = params.toList(par = o$par, spec = priors))
  
  if(hessian) {
    sigma.tmp = -solve(o$hessian)
    r$sigma =list(makePositiveDefinite(sigma.tmp[1:3,1:3]), 
                  makePositiveDefinite(sigma.tmp[-(1:3),-(1:3)]))
  }
  
  r
}