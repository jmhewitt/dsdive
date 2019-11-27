#' Convert vector of parameters on transformed scale to a list of parameters on model scale
#' 
#' The output also includes the log-Jacobian for the transformation.  Note that 
#' the transformation (and associated Jacobian) only links the unrestricted 
#' parameter space over which the Gibbs sampler is defined to the parameter 
#' space used by the model.  If prior distributions are defined for additional 
#' transformations of the model parameters, than the additional Jacobians must 
#' be accounted for within the \code{dsdive.prior} function.
#' 
#' @param par the vector of parameters to convert to a list of model parameters
#' @param spec List of parameters to specify prior distributions.  See 
#'   \code{dsdive.prior} for more details.
#'   
#' @import bisque
#' 
#' @example examples/params.convert.R
#' 
params.toList = function(par, spec) {
  
  # extract model parameters, and add jacobian for transformations
  res = list(
    # stage 1 downward trend, st. 2 directional persistence, st. 3 upward trend;
    # full matrix is used for backward-compatibility with an unconstrained, 
    # full-parameter model
    beta = matrix(c(bisque:::itx(x = par[1], link = 'logit', 
                        linkparams = list(range = c(0, spec$beta.absmax))), 
                    0, 0, par[2], 
                    bisque:::itx(x = par[3], link = 'logit', 
                        linkparams = list(range = c(-spec$beta.absmax, 0))),
                    0), nrow = 2),
    # transition rates, all stages
    lambda = exp(par[4:6]),
    # stage 1->2 transition parameters
    sub.tx = exp(par[7:8]),
    # stage 2->3 transition parameters
    surf.tx = exp(par[9:10]),
    # log-jacobian for transformations
    logJ = bisque:::jac.logit(x = par[1], log = TRUE, 
                              range = c(0, spec$beta.absmax)) + 
      bisque:::jac.logit(x = par[3], log = TRUE, 
                         range = c(-spec$beta.absmax, 0)) +
      sum(bisque:::jac.log(x = par[4:10], log = TRUE))
  )
  
  res 
}


#' Convert list of parameters to vector of parameters on transformed
#' 
#' @param par the list of parameters to convert to a vector of model parameters
#' @param spec List of parameters to specify prior distributions.  See 
#'   \code{dsdive.prior} for more details.
#'   
#' @import bisque
#' 
#' @example examples/params.convert.R
#'
params.toVec = function(par, spec) {
  c(bisque:::tx(x = par$beta[1,1], link = 'logit', 
                linkparams = list(range = c(0, spec$beta.absmax))),
    par$beta[2,2],
    bisque:::tx(x = par$beta[1,3], link = 'logit', 
                linkparams = list(range = c(-spec$beta.absmax, 0))),
    log(par$lambda), 
    log(par$sub.tx), 
    log(par$surf.tx)
  )
}