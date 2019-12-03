#' Evaluates prior distribution at parameter values
#' 
#' Note that the model parameters are converted to their unrestricted supports,
#' on which the prior distributions are defined.
#' 
#' @param par list specifying model parameters.  the list has entries: 
#'   \code{beta}, \code{lambda}, \code{sub.tx}, \code{surf.tx}
#' @param spec list specifying parameters for prior distributions.  the list has 
#'   entries:
#'   \describe{
#'     \item{beta.sd}{standard deviations for truncated normal prior 
#'       distributions associated with beta.}
#'     \item{beta.absmax}{maximum of absolute values for transformation out of 
#'       normal prior distributions associated with beta.}
#'     \item{lambda.sd}{standard deviations for centered normal prior 
#'       distributions associated with lambda.}
#'     \item{sub.tx.mean}{means for normal prior distributions associated with 
#'       stage 1->2 transition.}
#'     \item{sub.tx.sd}{standard devations for normal prior distributions 
#'       associated with stage 1->2 transition.}
#'     \item{surf.tx.mean}{means for normal prior distributions associated with 
#'       stage 2->3 transition.}
#'     \item{surf.tx.sd}{standard devations for normal prior distributions 
#'       associated with stage 2->3 transition.}
#'   }
#' @param log \code{TRUE} to return prior on log scale
#' 
#' @importFrom extraDistr dtnorm
#' @importFrom stats dnorm dlnorm
#' 
#' @export
#' 
#' @example examples/dsdive.prior.R
#' 
#'
dsdive.prior = function(par, spec, log = TRUE) {
  
  # soooo then, this step is not necessarily what we need since we *may* define 
  # prior distributions on a different support than the unrestricted space we 
  # use to move in the gibbs sampler.
  
  r = dtnorm(x = par$beta[1,1], sd = spec$beta.sd[1], log = TRUE, 
             a = 0, b = spec$beta.absmax) + 
    dnorm(x = par$beta[2,2], sd = spec$beta.sd[2], log = TRUE) +
    dtnorm(x = par$beta[1,3], sd = spec$beta.sd[3], log = TRUE, 
           a = -spec$beta.absmax, b = 0) + 
    sum(dlnorm(x = par$lambda, sdlog = spec$lambda.sd, log = TRUE)) +
    sum(dlnorm(x = par$sub.tx, meanlog = spec$sub.tx.mean, 
               sdlog = spec$sub.tx.sd, log = TRUE)) +
    sum(dlnorm(x = par$surf.tx, meanlog = spec$surf.tx.mean, 
               sdlog = spec$surf.tx.sd, log = TRUE))
  
  if(log) { r } else { exp(log) }
}