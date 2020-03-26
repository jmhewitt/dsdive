#' Likelihood for collection of partially observed dives
#' 
#' @param dsobs.list list of \code{dsobs} objects, which describe the 
#'   observation times and depths of a collection of dives
#' @param t.stages.list list of stage transition times for dives 
#'   observed in \code{dsobs.list}
#' @param P.raw list of continuous time probability transition matrices, and 
#'  components.
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
dsdive.obsld = function(dsobs.list, t.stages.list, P.raw, s0, sf) {
  
  # if a single dive is passed in, coerce it to list format
  if(inherits(dsobs.list, 'dsobs')) {
    dsobs.list = list(dsobs.list)
    t.stages.list = list(t.stages.list)
  }
  
  # compute likelihood for these stages
  s.range = s0:sf
  
  # number of depth bins
  n.bins = nrow(P.raw[[1]]$obstx.mat)
  
  # compute log-density for family of observations
  sum(sapply(1:length(dsobs.list), function(dive.id) {
    
    # extract observed dive components
    d = dsobs.list[[dive.id]]
    depths = d$depths
    times = d$times
    n = length(d$depths)
    
    # compute stage at each observation
    t.stages = t.stages.list[[dive.id]]
    stages = findInterval(times, t.stages) + 1
    
    ld = 0
    for(i in 1:(n-1)) {
      # extract start/end depth bins
      d0 = depths[i]
      df = depths[i+1]
      # extract start/end times
      t0 = times[i]
      tf = times[i+1]
      # extract start/end stages
      s0.step = stages[i]
      sf.step = stages[i+1]
      s.step = s0.step:sf.step
      
      # add likelihood contribution of observation
      if(any(s.step %in% s.range)) {
        
        # determine window of time spent in each stage
        dt.stages = sapply(s.step, function(s) {
          min(tf, t.stages[s], na.rm = TRUE) - max(t0, t.stages[s-1])
        })
        
        # look up or compute transition distribution through s0
        if(dt.stages[1] == P.raw[[s0.step]]$obstx.tstep) {
          u0 = P.raw[[s0.step]]$obstx.mat[d0,]
        } else {
          u0 = Matrix::expm(P.raw[[s0.step]]$A *  dt.stages[1])[d0,]
        }
        
        # diffuse transition distribution through other stages
        if(sf.step != s0.step) {
          for(s.ind in 2:length(s.step)) {
            s = s.step[s.ind]
            u0 = t(expm::expAtv(
              A = as.matrix(t(P.raw[[s]]$A)), 
              t = dt.stages[s.ind],
              v = as.numeric(u0)
            ))[[1]]
          }
        }
        
        # correct for numerical stability
        if(any(u0 < 0)) {
          # zap small negative values
          u0[(u0 < 0) & (u0 > -1e-10)] = 0
          # renormalize probability distribution
          u0 = u0 / sum(u0)
        }
        
        # add likelihood contribution
        ld = ld + log(u0[df])
      }
    }
    
    ld
  }))
  
}
