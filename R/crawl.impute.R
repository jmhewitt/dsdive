#' Use bridged sampling to impute a dive trajectory consistent with observations
#' 
#' The imputation will be done via crawl only, so the imputation only imputes 
#' likely locations in discrete space.  Dive stages will not be imputed.  For 
#' inference, the dive stages should be treated as missing data.
#'
#' 
#' @param s0 dive stage at which the trajectory should be started from
#' @param beta \eqn{2 x 3} matrix in which each column contains the diving 
#'  preference and directional persistence parameters for the DIVING, SUBMERGED, 
#'  and SURFACING dive stages.
#' @param lambda length 3 vector that specifies the transition rate, 
#'   respectively in the DIVING, SUBMERGED, and SURFACING stages.
#' @param sub.tx length 2 vector that specifies the first depth bin at which 
#'   transitions to the SUBMERGED stage can occur and the probability that such 
#'   a transition occurs at the next depth transition
#' @param surf.tx parameter that specifies the probability the trajectory will 
#'   transition to the SURFACING stage at the next depth transition
#' @param inflation.factor.lambda In order to facilitate bridged transitions, 
#'   the transition rate of the overall process must be inflated to allow the 
#'   possibility of self-transitions.  Self-transitions allow bridged paths to 
#'   dynamically modify the total number of transitions between observed values
#'   so that a valid path between observations is always possible.  The 
#'   \code{inflation.factor.lambda} parameter implicitly controls the number of 
#'   self-transitions that will occur.  Larger values will create more 
#'   self-transitions.
#' @param verbose If \code{TRUE}, then the sampler's progress will be printed 
#'   during sampling.
#' @param precompute.bridges If \code{TRUE}, then the bridged transition 
#'   matrices will be precomputed.  Enabling this option will increase the 
#'   memory overhead of the method, but will reduce its runtime.
#' 
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param depths record of depth bins the trajectory should visit
#' @param times times at which the depth bins should be visited
#' @param N number of trajectories to impute
#' @param depths.impute Either \code{'uniform'} or \code{'midpoint'}.  The 
#'   \code{crawl} model requires continuous rather than discrete depth values.
#'   Setting \code{depths.impute='uniform'} will fit the \code{crawl} model 
#'   using depth values that are uniformly sampled from each depth bin.  To 
#'   average over the uniform imputation, each of the \code{N} imputed 
#'   trajectories will fit the \code{crawl} model to a different set of 
#'   uniformly sampled depths.
#'   Setting \code{depths.impute='midpoint'} will fith the \code{crawl} model 
#'   using the midpoint of each depth bin.
#' @param tstep Time step to be used when imputing continuous trajectories to be 
#'   discretized
#' 
#' @example examples/crawl.impute.R
#' 
#' @importFrom crawl crwMLE crwSimulator crwPostIS
#' @import dplyr
#' 
#' @export
#'
crawl.impute = function(depth.bins, depths, times, N, 
                        depths.impute = 'uniform', tstep) {
  
  # precompute depth bin ranges
  depth.ranges = data.frame(min = depth.bins[,1] - depth.bins[,2],
                            max = depth.bins[,1] + depth.bins[,2])
  
  
  # impute continuous depths
  if(depths.impute == 'midpoint') {
    
    # use depth bin midpoint 
    df = data.frame(x = depth.bins[depths,1], time = times)
    
    # duplicate depth coordinates, for compatibility with crawl
    df$y = df$x
    
    # package in list
    depths.imputed = list(df)
    
  } else if(depths.impute == 'uniform') {
    
    # initialize storage
    depths.imputed = vector("list", length = N)
    
    # impute depths
    nt = length(depths)
    for(i in 1:N) {
      # use random depth within each observed bin
      df = data.frame( 
        x = sapply(depths, function(ind) {
          runif(n = 1, min = depth.ranges[ind,1], max = depth.ranges[ind,2])
          }),
        time = times
      )
    
      # duplicate depth coordinates, for compatibility with crawl
      df$y = df$x
      
      # package in list
      depths.imputed[[i]] = df
    }
    
  }
  
  
  # build times at which continuous depths should be imputed
  pred.times = seq(from = min(times), to = max(times), by = tstep)
  
  # number of paths each imputed depth record should generate
  paths.per.imputed = N / length(depths.imputed)
  
  
  # impute paths in continuous time and continuous depth
  sims = lapply(depths.imputed, function(df){
    
    # fit crawl model
    crawl.fit = crwMLE(data = df)
    
    # update variance to account for loss in data
    crawl.fit$Cmat = 2 * crawl.fit$Cmat
    crawl.fit$se = sqrt(diag(crawl.fit$Cmat))
    
    # create simulator
    sim.obj = crwSimulator(object.crwFit = crawl.fit, predTime = pred.times)
    
    # impute trajectories
    res = replicate(paths.per.imputed, list(
      crwPostIS(object.sim = sim.obj, fullPost = FALSE)
    ))
    
    # return result
    if(length(res)==1) {
      res[[1]]
    } else {
      res
    }
  })
  
  # remove a layer of depth from list structure
  if(N > 1) {
    if(length(sims)==1){ 
      sims = sims[[1]]
    }
  }
  
  # discretize paths
  sims.discretized = lapply(sims, function(s) {
      
    # discretize simulated depths to bins
    d = data.frame(depth.bin = findInterval(x = s$alpha.sim[,1],
                                            vec = depth.ranges$max) + 1,
                   t = s$time)
    
    # fill in missing transitions as necessary
    if(!all(unique(diff(d$depth.bin)) %in% c(-1,0,1))) {
      # identify indices where direct transitions were not observed
      inds = which(!(diff(d$depth.bin) %in% c(-1,0,1)))
      # splice in direct transitions; must be done in reverse order
      for(i in rev(inds)) {
        # linearly impute path
        path.missing = d$depth.bin[i]:d$depth.bin[i+1]
        path.missing = path.missing[-c(1,length(path.missing))]
        # linearly impute arrival times
        times.missing = seq(from = d$t[i], to = d$t[i+1], 
                            length.out = length(path.missing) + 2)
        times.missing = times.missing[-c(1,length(times.missing))]
        # splice in interpolation
        d = rbind(d[1:i,], cbind(depth.bin = path.missing, t =times.missing), 
                  d[(i+1):nrow(d),])
      }
    }
    
    # extract ctds
    res = as.list(d %>% 
                    mutate(ind = c(0, cumsum(abs(diff(depth.bin))))) %>% 
                    group_by(ind) %>%
                    summarise(depths = depth.bin[1], 
                              times = t[1]) %>% 
                    mutate(durations = c(diff(times), NA),
                           stages = -1) %>%
                    select(-ind))
    
    class(res) = 'dsdive'
    
    res
  })
  
  sims.discretized
}