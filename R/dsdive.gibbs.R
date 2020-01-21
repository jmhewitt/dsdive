checkForRemoteErrors = function (val, crash.fn) {
  count <- 0
  firstmsg <- NULL
  nodes = c()
  for(i in 1:length(val)) {
    v = val[[i]]
    if (inherits(v, c("try-error", 'character'))) {
      nodes = c(nodes, i)
      count <- count + 1
      if (count == 1)
        firstmsg <- v
    }
  }
  if(count > 0 ) {
    print(firstmsg)
    crash.fn(nodes, firstmsg)
  }
  if (count == 1)
    stop("one node produced an error: ", firstmsg)
  else if (count > 1)
    stop(count, " nodes produced errors; first error: ",
         firstmsg)
  val
}

#' Likelihood for completely observed dive trajectories
#'
#' @param depths Depth bin indices in which the trajectory was observed
#' @param durations Amount of time spent in each depth bin
#' @param times Times at which each depth bin was entered
#' @param stages Record of dive stages at each depth bin
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
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param t0.dive Time at which dive started
#' @param d0.last If the depth bin that proceeded the first depth bin in 
#'   \code{depths}.  If the trajectory to be analyzed was started at the 
#'   surface, then set \code{c0.last=NULL}.
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#' 
#' @example examples/ld.R
#' 
#' @importFrom snow checkCluster sendCall recvResult
#' 
#' 
dsdive.gibbs = function(
  dives.obs, cl, impute.init, impute.gibbs, init, verbose = FALSE, maxit, 
  checkpoint.fn, checkpoint.interval, T1.prior, T2.prior, pi1.prior, pi2.prior, 
  lambda1.prior, lambda2.prior, lambda3.prior, crash.fn) {
  
  #
  # initialize cluster
  #
  
  if(verbose) {
    message('Initializing cluster')
  }
  
  # validate cluster
  checkCluster(cl)
  
  # set up basic load-balancing scheme to partition data across cluster
  batch.size = floor(length(dives.obs) / length(cl))
  remainder = length(dives.obs) - batch.size * length(cl)
  
  # initialize cluster
  nodes.used = NULL
  start = 1
  for(i in seq(along = cl)) {
    
    # determine which dives go to node i
    end = start + batch.size + ifelse(i <= remainder, 1, 0) - 1
    
    # initialize node i
    if(end >= start) {
      
      inds = start:end
      nodes.used = c(nodes.used, i)
      
      if(verbose) {
        message(paste('    Initializing node', i, sep = ' '))
      }
      
      local({
        env <- as.environment(1) ## .GlobalEnv
        gets <- function(n, v) { assign(n, v, envir = env); NULL }
        sendCall(con = cl[[i]], fun = function(pkg, assign.fn, impute.init, 
                                               theta) {
          
          # ensure node has loaded tools
          library(dsdive)
          
          # impute initial trajectories
          imputed.local = lapply(pkg$dives.obs, function(d) {
            
            # initialize stage duration
            t.stages = seq(from = d$dive$times[1], 
                           to = d$dive$times[length(d$dive$times)],
                           length.out = 4)[2:3]
            
            # impute trajectory
            r = impute.init(depth.bins = d$depth.bins, depths = d$dive$depths,
                           times = d$dive$times, beta = theta$beta,
                           lambda = theta$lambda, t.stages = t.stages)
            
            r$t.stages = t.stages
            r
          })

          # save data to node
          assign.fn('impute.gibbs', pkg$impute.gibbs)
          assign.fn('T1.prior', pkg$T1.prior)
          assign.fn('T2.prior', pkg$T2.prior)
          assign.fn('obs.local', pkg$dives.obs)
          assign.fn('imputed.local', imputed.local)
          assign.fn('assign.fn', assign.fn)
          
          # return imputed trajectories
          list(imputed.local)
          
        }, args = list(pkg = list(dives.obs = dives.obs[inds], 
                                  T1.prior = T1.prior, T2.prior = T2.prior, 
                                  impute.gibbs = impute.gibbs), 
                       assign.fn = gets, impute.init = impute.init, 
                       theta = init))
      })
      
    }
    start = end + 1
  }
  
  # collect initial trajectory imputations
  dives.imputed = checkForRemoteErrors(sapply(cl[nodes.used], recvResult), 
                                       crash.fn)
  dives.imputed = do.call(c, dives.imputed)
  
  if(verbose) {
    message('Cluster initialized')
  }
  

  #
  # initialize sampler state and output
  #
  
  theta = init
  
  trace = matrix(NA, nrow = maxit, ncol = 5)
  colnames(trace) = c('pi1', 'pi2', 'lambda1', 'lambda2', 'lambda3')
  imputed.trace = vector('list', maxit)
  
  tick.checkpoint = proc.time()[3]
  
  
  #
  # gibbs sample
  #
  
  for(it in 1:maxit) {
    
    if(verbose) {
      message(paste('Iteration:', it, sep = ' '))
      tick = proc.time()[3]
    }
    
    
    #
    # impute trajectories
    #
    
    for(i in nodes.used) {
      local({
        env <- as.environment(1) ## .GlobalEnv
        gets <- function(n, v) { assign(n, v, envir = env); NULL }
        sendCall(con = cl[[i]], fun = function(theta) {
          for(j in 1:length(imputed.local)) {
            
            d = obs.local[[j]]
            r = imputed.local[[j]]
            t.stages = r$t.stages
            
            # impute trajectory
            r = impute.gibbs(depth.bins = d$depth.bins, depths = d$dive$depths, 
                             times = d$dive$times, beta = theta$beta, 
                             lambda = theta$lambda, imputed.cond = r)
            
            # sample stage transition times and updated stage vector
            stages.cond = dsdive.sample.stages(
              depths = r$depths, durations = r$durations, times = r$times,
              t.stages = t.stages, beta = theta$beta, lambda = theta$lambda,
              depth.bins = d$depth.bins, T1.prior = T1.prior, 
              T2.prior = T2.prior
            )
            
            r$stages = stages.cond$stages
            r$t.stages = stages.cond$t.stages
            r$suffstat = dsdive.suffstat(depths = r$depths, 
                                         durations = r$durations, 
                                         stages = r$stages, 
                                         depth.bins = d$depth.bins)
            
            imputed.local[[j]] = r
          }
          
          # update imputed dives on node
          assign.fn('imputed.local', imputed.local)
          
          # return imputed trajectories
          imputed.local
          
        }, args = list(theta = theta))
      })
    }
    
    # collect trajectory imputations
    dives.imputed = checkForRemoteErrors(sapply(cl[nodes.used], recvResult), 
                                         crash.fn)
    dives.imputed = do.call(c, dives.imputed)
    
    #
    # sample model parameters from full conditional posterior densities
    #
    
    dat.all = sapply(dives.imputed, function(r) {
      unlist(r$suffstat)
    })
    
    suffstats.all = rowSums(dat.all)
    
    theta$beta = c(
      rbeta(n = 1, 
            shape1 = pi1.prior[1] + suffstats.all['n.down1'], 
            shape2 = pi1.prior[2] + suffstats.all['n.beta1'] - 
              suffstats.all['n.down1']),
      rbeta(n = 1, 
            shape1 = pi2.prior[1] + suffstats.all['n.down2'], 
            shape2 = pi2.prior[2] + suffstats.all['n.beta2'] - 
              suffstats.all['n.down2'])
    )
    
    theta$lambda = c(
      rgamma(n = 1, 
             shape = lambda1.prior[1] + suffstats.all['n.lambda1'], 
             rate = lambda1.prior[2] + suffstats.all['d.lambda1']),
      rgamma(n = 1, 
             shape = lambda2.prior[1] + suffstats.all['n.lambda2'], 
             rate = lambda2.prior[2] + suffstats.all['d.lambda2']),
      rgamma(n = 1, 
             shape = lambda3.prior[1] + suffstats.all['n.lambda3'], 
             rate = lambda3.prior[2] + suffstats.all['d.lambda3'])
    )
    
    if(verbose) {
      print(theta)
    }
    
    
    #
    # save trace
    #
    
    # save trace 
    trace[it,1:5] = unlist(theta[1:2])
    imputed.trace[[it]] = dives.imputed
    
    tock = proc.time()[3]
    
    if(verbose) {
      message(paste('   ', round(tock - tick, 1), 'seconds', sep = ' '))
    }
    
    
    #
    # checkpoint behaviors
    #
    
    # dump state
    if(tock - tick.checkpoint > checkpoint.interval) {
      if(verbose) {
        message('--- Checkpoint ---')
      }
      checkpoint.fn(trace = trace[1:it,], imputed.trace = imputed.trace[1:it])
      tick.checkpoint = proc.time()[3]
    }
    
  }
  
  list(
    theta = trace,
    dives = imputed.trace
  )
}