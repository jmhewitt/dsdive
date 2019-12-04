#' Define generic method for evaluating dive model likelihoods under different
#' computing environments
#'
#' 
dsdive_ld = function(cfg, params) {
  UseMethod("dsdive_ld")
}