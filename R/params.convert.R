#' Convert vector of parameters on transformed scale to a list of parameters on model scale
#' 
#' The output also includes the log-Jacobian for the transformation
#' 
#' @param par the vector of parameters to convert to a list of model parameters
#' @param sub.tx1 the fixed value of sub.tx1 to incorporate into the parameter 
#'   list
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
params.toList = function(par, sub.tx1) {
  
  # extract model parameters, and add jacobian for transformations
  res = list(
    beta = matrix(par[1:6], nrow = 2),
    lambda = exp(par[7:9]),
    sub.tx = c(sub.tx1, plogis(par[10])),
    surf.tx = par[11:12],
    logJ = jac.logit(x = par[10], log = TRUE) + 
      sum(jac.log(x = par[7:9], log = TRUE))
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
  c(par$beta, log(par$lambda), qlogis(par$sub.tx[2]), par$surf.tx)
}

