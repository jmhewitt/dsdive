#' Define generic method for evaluating dive model likelihoods under different
#' computing environments
#'
#' All log-densities are evaluted with respect to the parameters on their 
#' natural scales, vs. log-densities with respect to the transformed scales.
#' 
dsdive_ld = function(cfg, params) {
  UseMethod("dsdive_ld")
}