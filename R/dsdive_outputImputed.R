#' Define generic method for backing up gibbs sampler state under different
#' computing environments
#'
#' 
dsdive_outputImputed = function(cfg, output, save.time, file) {
  UseMethod("dsdive_outputImputed")
}