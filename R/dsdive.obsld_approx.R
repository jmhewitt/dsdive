#' Likelihood for collection of partially observed dives
#' 
#' @param dsobs.list list of \code{dsobs} objects, which describe the 
#'   observation times and depths of a collection of dives
#' @param t.stages.list list of stage transition times for dives 
#'   observed in \code{dsobs.list}
#' @param P.interpolator list of functions to approximate transition 
#'   probabilities in a given stage, given model parameters, and a timestep
#' @param beta vector with directional preference parameters for descent and  
#'   ascent stages.
#' @param lambda vector with dive rate parameters for the descent, bottom, and 
#'   ascent stages.
#' @param s0 likelihood should include contributions for stages greater or 
#'   equal to \code{s0}
#' @param sf likelihood should include contributions for stages less than or 
#'   equal to \code{sf}
#'
#' @importFrom Matrix sparseVector expm
#' @importFrom expm expAtv
#' 
#' @example examples/dsdive.obsld.R
#' 
#' @export
#'
dsdive.obsld_approx = function(dsobs.list, t.stages.list, s0, sf, beta, lambda, 
                               P.interpolator) {
  
  # if a single dive is passed in, coerce it to list format
  if(inherits(dsobs.list, 'dsobs')) {
    dsobs.list = list(dsobs.list)
    t.stages.list = list(t.stages.list)
  }
  
  # compute likelihood for these stages
  s.range = s0:sf
  
  # compute log-density for family of observations
  sum(sapply(1:length(dsobs.list), function(dive.id) {
    
    # extract observed dive components
    d = dsobs.list[[dive.id]]
    depths = d$depths
    times = d$times
    n = length(d$depths)
    
    # extract exact stage transition times
    t.stages = t.stages.list[[dive.id]]
    
    # round stage transition times to closest observation time to avoid dense 
    # matrix operations required by arbitrary stage transition times
    t.stages = times[apply(abs(outer(times, t.stages, `-`)), 2, which.min)]
    
    # compute stage at each observation
    stages = findInterval(times, t.stages) + 1
    
    ld = 0
    for(i in 1:(n-1)) {
      # extract start/end depth bins
      d0 = depths[i]
      df = depths[i+1]
      # extract duration
      dt = times[i+1] - times[i]
      # extract start stage
      s0.step = stages[i]
      
      # approximate likelihood contribution of observation
      if(s0.step %in% s.range) {
        # extract stage parameters
        if(s0.step == 1) {
          beta.step = beta[1]
          lambda.step = lambda[1]
        } else if(s0.step == 2) {
          beta.step = .5
          lambda.step = lambda[2]
        } else if(s0.step == 3) {
          beta.step = beta[2]
          lambda.step = lambda[3]
        }
        # approximate transition probability
        log_p = P.interpolator[[s0.step]](
          beta = beta.step, lambda = lambda.step, tstep = dt, 
          i = d0, j = df, log = TRUE)
        
        # aggregate likelihood contribution
        ld = ld + log_p
      }
    }
    
    ld
  }))
  
}
