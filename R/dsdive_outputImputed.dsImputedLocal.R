#' Define generic method for backing up gibbs sampler state under different
#' computing environments
#'
#' 
dsdive_outputImputed.dsImputedLocal = function(cfg, output, save.time, 
                                               file) {
  
  # save local copies of trajectories
  lapply(cfg$ids, function(id, file) {
    
    # extract local config
    cfg.local = get(x = id, envir = globalenv())
    
    # build and save to local file
    file.local = paste(gsub(pattern = '\\.RData', replacement = '', x = file),
                       '_imputed_trace_', id, '.RData', sep = '')
    save(cfg.local, file = file.local)
    
  }, file = file)
  
  # return output unmodified
  output
}