#' Compute probability transition matrix for partially observed CTMC
#' 
#' Given model parameters, this function will compute the probability transition 
#' matrix for a Continuous time Markov chain (CTMC) that is observed once, and 
#' then again at a time of \code{tstep} units of time later.  If 
#' \code{include.raw==TRUE}, then raw components of the probability transition 
#' matrix will be returned so that the transition matrix can be computed for 
#' arbitrary timesteps.
#' 
#' @param pi.designs list of design matrices for the diving preference 
#'   parameters by stage.
#' @param lambda.designs list of design matrices for the diving speed
#'   parameters by stage.
#' @param beta1 coefficient values for linear model for logit(pi^{(1)}).
#' @param beta2 coefficient values for linear model for logit(pi^{(3)}).
#' @param alpha1 coefficient values for linear model for log(lambda^{(1)}).
#' @param alpha2 coefficient values for linear model for log(lambda^{(2)}).
#' @param alpha3 coefficient values for linear model for log(lambda^{(3)}).
#' @param s0 stage for which the transition matrix should be computed
#' @param ind transition matrix is unique for each dive.  The parameter 
#'   \code{ind} specifies that the transition matrix should be computed for the 
#'   dive associated with the covariates in row \code{ind} in the design 
#'   matrices (see \code{pi.designs}, \code{lambda.designs}).
#' @param tstep Time between observations of the CTMC
#' @param include.raw \code{TRUE} to include raw components of the probability
#'   transition matrix so that the matrix can be computed for arbitrary 
#'   timesteps
#' @param delta If \code{delta>0}, then the observation matrix and raw 
#'   components computed will be for a transition matrix whose generator is 
#'   perturbed to allow much faster computation.
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' 
#' @export
#' 
#' @importFrom stats plogis
#' 
#' @example examples/dsdive.obstxmat.cov.R
#' 
dsdive.obstxmat.cov = function(pi.designs, lambda.designs, beta1, beta2, 
                               alpha1, alpha2, alpha3, s0, ind, 
                               tstep, include.raw, depth.bins, delta) {
  
  #
  # compute dive-specific parameters on transformed scale for stage s0
  #
  
  pi = numeric(2)
  lambda = numeric(3)
  
  if(s0==1) {
    pi[1] = pi.designs[[1]][ind,, drop = FALSE] %*% beta1[]
    lambda[1] = lambda.designs[[1]][ind,, drop = FALSE] %*% alpha1[]
  } else if(s0==2) {
    lambda[2] = lambda.designs[[2]][ind,, drop = FALSE] %*% alpha2[]
  }else if(s0==3) {
    pi[2] = pi.designs[[2]][ind,, drop = FALSE] %*% beta2[]
    lambda[3] = lambda.designs[[3]][ind,, drop = FALSE] %*% alpha3[]
  }
  
  # back-transform parameters
  pi = plogis(pi)
  lambda = exp(lambda)
  
  # build requested transition matrix
  dsdive.obstx.matrix(depth.bins = depth.bins, beta = pi, lambda = lambda, 
                      s0 = s0, tstep = tstep, include.raw = include.raw, 
                      delta = delta)
}