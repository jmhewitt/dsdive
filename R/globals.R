#' @useDynLib dsdive, .registration = TRUE
globalVariables(c('ind','depth.bin', 't.start', 't.end', 'depth.min', 
                  'depth.max', 'stages', 'depth.mid', 'local_n', 'local_ids', 
                  'beta1', 'beta2', 'alpha1', 'alpha2', 'alpha3', 'pi.designs', 
                  'lambda.designs', 'tstep', 'depth.bins', 'delta', 
                  'dsobs.list', 't.stages', 'myinfo'))