#' Sample stage transition times from full conditional posterior
#' 
#' @param prior.t list of two functions that will allow the prior density for 
#'   stage transition times to be called
#' 
#' @example examples/dsdive.sample.stages.R
#' 
#' @export
#' 
dsdive.sample.stages = function(depths, times, t.stages, beta, lambda, 
                                depth.bins, T1.logprior, T2.logprior, t1.sd, 
                                t2.sd) {
  
  # compute stages, durations, and log-density for input
  
  augmented = dsdive.augment.trajectory(depths = depths, times = times, 
                                        t.stages = t.stages)
  
  ld = dsdive.ld.fixedstages(depths = augmented$depths, 
                             durations = augmented$durations, 
                             times = augmented$times, 
                             stages = augmented$stages, 
                             beta = beta, lambda = lambda, 
                             depth.bins = depth.bins)
    
  
  #
  # sample stage 1->2 transition time, conditional on stage 3 transition time
  #
  
  # get start time for dive
  t0.dive = times[1]
  
  # propose new transition time
  T1.prop = proposal.logit(x0 = t.stages[1], sd = t1.sd, 
                           a = t0.dive, b = t.stages[2])
  
  # compute stages, durations, and log-density for proposal
  
  augmented.prop = dsdive.augment.trajectory(depths = depths, times = times, 
                                             t.stages = c(T1.prop$x, 
                                                          t.stages[2]))
  
  ld.prop = dsdive.ld.fixedstages(depths = augmented.prop$depths, 
                                  durations = augmented.prop$durations, 
                                  times = augmented.prop$times, 
                                  stages = augmented.prop$stages, 
                                  beta = beta, lambda = lambda, 
                                  depth.bins = depth.bins)
  
  # accept/reject
  lR = ld.prop + T1.logprior((T1.prop$x - t0.dive)/60) - ld - 
    T1.logprior((t.stages[1] - t0.dive)/60) + T1.prop$lR
  a1 = log(runif(1)) < lR
  if(a1) {
    t.stages[1] = T1.prop$x
    ld = ld.prop
    augmented = augmented.prop
  }
  
  
  #
  # sample stage 2->3 transition time, conditional on stage 2 transition time
  #
  
  # get end time for dive
  tf.dive = times[length(times)]
  
  # propose new transition time
  T2.prop = proposal.logit(x0 = t.stages[2], sd = t2.sd, 
                           a = t.stages[1], b = tf.dive)
  
  # compute stages, durations, and log-density for proposal
  
  augmented.prop = dsdive.augment.trajectory(depths = depths, times = times, 
                                             t.stages = c(t.stages[1], 
                                                          T2.prop$x))
  
  ld.prop = dsdive.ld.fixedstages(depths = augmented.prop$depths, 
                                  durations = augmented.prop$durations, 
                                  times = augmented.prop$times, 
                                  stages = augmented.prop$stages, 
                                  beta = beta, lambda = lambda, 
                                  depth.bins = depth.bins)
  
  # accept/reject
  lR = ld.prop + T2.logprior((T2.prop$x - t.stages[1])/60) - ld - 
    T2.logprior((t.stages[2] - t.stages[1])/60) + T2.prop$lR
  a2 = log(runif(1)) < lR
  if(a2) {
    t.stages[2] = T2.prop$x
    ld = ld.prop
    augmented = augmented.prop
  }
  
  simplified = dsdive.simplify.trajectory(depths = augmented$depths, 
                                          times = augmented$times, 
                                          stages = augmented$stages)
  
  # package results
  list(
    dive = simplified,
    ld = ld, 
    t.stages = t.stages,
    accepted = c(a1, a2)
  )
}