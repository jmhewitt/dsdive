#' Compute probability transition matrix for partially observed CTMC
#' 
#' Given model parameters, this function will compute the probability transition 
#' matrix for a Continuous time Markov chain (CTMC) that is observed once, and 
#' then again at a time of \code{tstep} units of time later.  If 
#' \code{include.raw==TRUE}, then raw components of the probability transition 
#' matrix will be returned so that the transition matrix can be computed for 
#' arbitrary timesteps.
#' 
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param beta the current depth bin transition model parameters
#' @param lambda the current depth bin transition rate model parameters
#' @param s0 the stage for which to compute the transition matrix
#' @param tstep Time between observations of the CTMC
#' @param include.raw \code{TRUE} to include raw components of the probability
#'   transition matrix so that the matrix can be computed for arbitrary 
#'   timesteps
#'   
#' @example examples/dsdive.obstx.matrix.R
#' 
#' @importFrom Matrix expm
#' 
#' @export
#' 
dsdive.obstx.matrix = function(depth.bins, beta, lambda, s0, tstep, 
                               include.raw = FALSE) {
  
  # build uniformized generator matrix
  rate.unif = max(lambda[s0] / (2*depth.bins[,2]))
  A = dsdive.generator.matrix.uniformized(
    depth.bins = depth.bins, beta = beta, lambda = lambda, s0 = s0, 
    rate.uniformized = rate.unif)
  
  # compute matrix exponential
  obstx.mat = Matrix::expm(A * tstep)

  if(include.raw) {
    list(
      obstx.mat = obstx.mat,
      obstx.tstep = tstep,
      A = A
    )
  } else {
    obstx.mat
  }
}