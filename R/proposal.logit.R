#' Logit-based Metropolis-Hastings Random-Walk proposals
#' 
#' Uses a logit transform on x0 to propose a new parameter value, for use in 
#' MH RW Gibbs steps.  The function also returns the log-ratio of the proposal.
#' 
#' @example examples/proposal.logit.R
#' 
#' @param x0 initial parameter value
#' @param sd standard deviation of proposal, on logit scale
#' @param a minimum allowable value for the proposal
#' 
#' 
proposal.logit = function(x0, sd, a = 0, b = 1) {
  
  # get width of support for x
  w = b-a
  
  # transform x0 to logit scale
  u0 = qlogis((x0-a)/w)
  
  # propose new value of x on logit scale
  u = u0 + rnorm(n = 1, sd = sd)
  
  # back-transform proposal
  x = a + w * plogis(u)
  
  # compute log-ratio for proposal
  lR = log(x-a) + log(b-x) - log(x0-a) - log(b-x0)
  
  list(x = x, lR = lR)
}