#' Define generic method for backing up gibbs sampler state under different
#' computing environments
#'
#' 
#' 
dsdive_outputImputed.dsCompleteDistributed = function(cfg, output, save.time, 
                                                      file) {
  # no trajectories were imputed, so no backup actions required
  output
}