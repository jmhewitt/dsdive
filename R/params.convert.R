#' Convert vector of parameters on transformed scale to a list of parameters on model scale
#' 
#' The output also includes the log-Jacobian for the transformation
#' 
#' @param par the vector of parameters to convert to a list of model parameters
#'   
#' @importFrom bisque jac.logit jac.log
#' 
#' @examples 
#' 
#' data('dive.sim')
#' attach(dive.sim)
#' attach(params)
#' 
#' params.toList(par = runif(12), sub.tx1 = 1)
#' 
#' detach(params)
#' detach(dive.sim)
#' 
params.toList = function(par) {
  
  # extract model parameters, and add jacobian for transformations
  res = list(
    # stage 1 downward trend, st. 2 directional persistence, st. 3 upward trend;
    # full matrix is used for backward-compatibility with an unconstrained, 
    # full-parameter model
    beta = matrix(c(exp(par[1]), 0,
                    0, par[2], 
                    -exp(par[3]), 0), nrow = 2),
    # transition rates, all stages
    lambda = exp(par[4:6]),
    # stage 1->2 transition parameters
    sub.tx = par[7:8],
    # stage 2->3 transition parameters
    surf.tx = par[9:10],
    # log-jacobian for transformations
    logJ = sum(jac.log(x = par[c(1,3,4:6)], log = TRUE))
  )
  
  res 
}


#' Convert list of parameters to vector of parameters on transformed
#' 
#' @param par the list of parameters to convert to a vector of model parameters
#' 
#' @examples 
#' 
#' data('dive.sim')
#' attach(dive.sim)
#' attach(params)
#' 
#' params.toVec(par = params)
#' 
#' detach(params)
#' detach(dive.sim)
#'
params.toVec = function(par) {
  c(log(par$beta[1,1]), par$beta[2,2], log(-par$beta[1,3]), log(par$lambda),
    par$sub.tx, par$surf.tx)
}