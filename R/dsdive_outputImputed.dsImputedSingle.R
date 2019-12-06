#' Define generic method for backing up gibbs sampler state under different
#' computing environments
#'
#' 
#' 
dsdive_outputImputed.dsImputedSingle = function(cfg, output, save.time, file) {
  # extract and return imputed paths
  output$trace.imputed = cfg$trace.imputed
  output
}