#' Define generic method for backing up gibbs sampler state under different
#' computing environments
#'
#' @importFrom snow clusterApply
#' 
dsdive_outputImputed.dsImputedBatchDistributed = function(cfg, output, 
                                                          save.time, file) {

  # save local copies of trajectories
  clusterApply(cl = cfg$cl, x = cfg$cluster.ids, function(id, file) {
    
    # extract local config
    cfg.local = get(x = id, envir = globalenv())
    
    if(!is.null(cfg.local)) {
      file.batch = paste(gsub(pattern = '\\.RData', replacement = '', x = file),
                         '_batch_', id, '.RData', sep = '')
      dsdive_outputImputed(cfg = cfg.local, output = output, 
                           save.time = save.time, file = file.batch)
    }
    
  }, file = file)
  
  # return output unmodified
  output
}