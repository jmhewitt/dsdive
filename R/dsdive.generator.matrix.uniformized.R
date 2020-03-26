#' Compute infinitesimal generator matrix for a homogeneous, uniformized CTMC
#' 
#' Compute infinitesimal generator matrix for a uniformized Continuous time
#' Markov chain (CTMC).  See Rao and Teh (2013) for an applied definition and 
#' use of infinitesimal generator matrices.
#' 
#' @references Rao, Vinayak, and Yee Whye Teh. "Fast MCMC sampling for Markov 
#'   jump processes and extensions." The Journal of Machine Learning Research 
#'   14.1 (2013): 3295-3320.
#'   
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param beta Directional preference model parameters.  See 
#'   \code{dsdive.tx.params} for more details.
#' @param lambda Diving rate model parameters.  See 
#'   \code{dsdive.tx.params} for more details.
#' @param s0 dive stage for which matrix should be computed
#' @param rate.uniformized uniformization rate, for standardizing transition
#'   rates between states
#'   
#' @example examples/dsdive.generator.matrix.uniformized.R
#' 
#' @importFrom Matrix diag rowSums
#' 
#' @export
#' 
dsdive.generator.matrix.uniformized = function(depth.bins, beta, lambda, s0,
                                               rate.uniformized) {
  
  # start by building uniformized transition matrix
  A = dsdive.tx.matrix.uniformized(
    depth.bins = depth.bins, beta = beta, lambda = lambda, s0 = s0, 
    rate.uniformized = rate.uniformized)
  
  
  # rescale matrix and diagonal entries
  A = A * rate.uniformized
  diag(A) = diag(A) - rate.uniformized
  
  A
}