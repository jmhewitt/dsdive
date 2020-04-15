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
#' @param A intended to be used to pass in pre-allocated, shared-memory storage 
#'   for the uniformized transition matrix.  Shared memory is managed by the 
#'   \code{bigmemory} package.  If \code{A=NULL}, then a standard R matrix will
#'   be initialized and returned instead.
#' @param obstx.mat intended to be used to pass in pre-allocated, shared-memory 
#'   storage for the probability transition matrix.  Shared memory is managed 
#'   by the \code{bigmemory} package.  If \code{obstx.mat=NULL}, then a 
#'   standard R matrix will be initialized and returned instead.
#'   
#' @example examples/dsdive.obstx.matrix.R
#' 
#' @importFrom Matrix expm
#' @import bigmemory
#' 
#' @useDynLib dsdive, .registration = TRUE
#' 
#' @export
#' 
dsdive.obstx.matrix = function(depth.bins, beta, lambda, s0, tstep, 
                               include.raw = FALSE, A = NULL, 
                               obstx.mat = NULL) {
  
  # build uniformized generator matrix
  rate.unif = max(lambda[s0] / (2*depth.bins[,2]))
  A = dsdive.generator.matrix.uniformized(
    depth.bins = depth.bins, beta = beta, lambda = lambda, s0 = s0, 
    rate.uniformized = rate.unif, A = A)
  
  # compute matrix exponential
  res = expm_cpp(A = as.matrix(A[]), delta = 1e-10, t = tstep)
  if(is.null(obstx.mat)) {
    # obstx.mat = Matrix::expm(A[] * tstep)
    obstx.mat = res$expm
  } else {
    # obstx.mat[] = Matrix::expm(A[] * tstep)
    obstx.mat[] = res$expm
  }
  
  if(include.raw) {
    list(
      obstx.mat = obstx.mat,
      obstx.tstep = tstep,
      A = A,
      vectors = res$vectors,
      values = res$values,
      d = res$d,
      dInv = res$dInv
    )
  } else {
    obstx.mat
  }
}