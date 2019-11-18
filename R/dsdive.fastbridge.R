#' Impute a family of dive trajectories that meet start and end times
#'
#' Sample a dive trajectory such that its start and end locations and times are 
#' fixed.  Sampling uses a thinned Poisson process and bridged Markov chain 
#' transition probabilities.  Sampling is not from the model's theoretical 
#' imputation distribution, but an approximation thereof.  The trajectory's
#' log-density is returned with the output.
#' 
#' The sampling method is designed to sample many trajectories simultaneously, 
#' so has an extra level of approximation in the proposal distributions.
#' 
#' @param M the number of trajectories to sample
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
#' @param d0 the depth bin at which the trajectory begins
#' @param d0.last the depth bin from which the trajectory entered \code{d0}.
#'   Set \code{d0=NULL} if the trajectory is beginning at the surface depth bin.
#' @param df the depth bin at which the trajectory should end
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
#' @param t0 time at which the trajectory should start
#' @param tf time at which the trajectory should end
#' @param s0 dive stage in which the trajectory begins
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
#' @param lambda.max Arrival rate for the parent Poisson process that will
#'   be thinned.  \code{lambda.max} will be scaled by 
#'   \code{inflation.factor.lambda}, and if \code{lambda.max==NULL} then the 
#'   method will compute this on its own.
#' @param t0.dive Time at which dive started
#' @param trajectory.conditional If not \code{NULL}, then 
#'   \code{trajectory.conditional} must be the dive information for a 
#'   completely observed \code{dsdive} object.  The first entries in the 
#'   initialization vectors \code{d0.last} and \code{s0} must be 
#'   associated with the trajectory observed in \code{trajectory.conditional}.
#'   Providing a non \code{NULL} value for \code{trajectory.conditional} will 
#'   cause \code{dsdive.fastbridge} to simulate one fewer trajectories, as 
#'   the value of \code{trajectory.conditional} will be returned.  The 
#'   importance of this function argument is that \code{dsdive.fastbridge}
#'   will evaluate the proposal density for \code{trajectory.conditional}.
#' 
#' @example examples/dsdive.fastbridge.R
#' 
#' @importFrom extraDistr rtpois dtpois
#' @import Matrix
#' 
#' @export
#' 
#'
dsdive.fastbridge = function(M, depth.bins, d0, d0.last, df, beta, lambda, 
                             sub.tx, surf.tx, t0, tf, s0, 
                             inflation.factor.lambda = 1.1, verbose = FALSE,
                             precompute.bridges = TRUE, lambda.max = NULL,
                             t0.dive, trajectory.conditional = NULL) {
  
  #
  # build basic simulation parameters
  #
  
  # determine the actual number of trajectories that must be sampled
  cond.sim = !is.null(trajectory.conditional)
  if(cond.sim) {
    M.sim = M - 1
  } else { 
    M.sim = M
  }
  
  # initialize log-proposal density for samples
  ld = numeric(M)
  
  # transition time window
  T.win = tf - t0
  
  # minimum number of transitions required for bridging
  min.tx = abs(df-d0)
  
  # convert d0.last to vector if not already provided
  if(length(d0.last) != M) {
    d0.last = rep(d0.last[1], M)
  }
  
  # convert s0 to vector if not already provided
  if(length(s0) != M) {
    s0 = rep(s0[1], M)
  }
  
  
  #
  # sample homogeneous poisson process for thinning 
  #
  
  # get rate for thick poisson process
  if(is.null(lambda.max)) {
    lambda.max = max(outer(lambda, 2 * depth.bins[,2], '/'))
  }
  lambda.thick = inflation.factor.lambda * lambda.max
    
  
  # sample number of unthinned arrivals
  lambda.tmp = T.win * lambda.thick
  N = rtpois(n = M.sim, lambda = lambda.tmp, a = min.tx)
  
  # preprend the number of transitions for trajectory.conditional
  if(cond.sim) {
    N = c(length(trajectory.conditional$depths)-1, N)
  }
  
  # quick fix to resample infinite draws
  N.inf = is.infinite(N)
  if(any(N.inf)) {  
    N[N.inf] = sample(x = N[-N.inf], size = sum(N.inf), replace = TRUE)
  } 
  
  # get most number of arrivals required during sampling
  N.max = max(N)
  
  
  if(verbose) {
    message(paste('Simulating total of', sum(N), 'potential transitions', 
                  sep =' '))
  }
  
  # update log-proposal density for sample
  ld = ld +
    # log-density for number of arrivals
    # note: a = min.tx - 1 accounts for bug in dtpois support
    dtpois(x = N, lambda = lambda.tmp, a = min.tx-1, log = TRUE) +
    # log-density for arrival times (as order statistics of uniform sample)
    lfactorial(N) - N * log(T.win)
    
  
  #
  # bridge sample trajectories
  #
    
    if(verbose) {
      message('Computing transition matrices at potential transition times')
    }
    
    # build complete transition matrices for each arrival time
    min.depth = max(1, d0 - N.max)
    max.depth = min(nrow(depth.bins), d0 + N.max)
    
    # compute the base transition matrix
    #
    # note: in ideal cases, the transition matrix only has a weak dependence on 
    # t0, so we accelerate the sampling process by basing all proposals off of a
    # single transition matrix
    tx.mat.raw = dsdive.tx.matrix(t0 = t0, depth.bins = depth.bins, beta = beta, 
                                  lambda = lambda, sub.tx = sub.tx, 
                                  surf.tx = surf.tx,
                                  inflation.factor.lambda = 1,
                                  min.depth = min.depth, max.depth = max.depth,
                                  t0.dive = t0.dive, lambda.max = lambda.thick)
    tx.mat = tx.mat.raw$m
    out.inds.lookup = tx.mat.raw$out.inds
    
    # add a "null" depth bin to allow trajectory initialization
    n = nrow(depth.bins) + 1
    
    # get ending indices for trajectory (any previous depth, any dive stage)
    end.inds = c()
    for(s in min(s0):3) {
      for(dd in c(-1,0,1)) {
        if(((df + dd) > 0) & ((df + dd) < n)) {
          end.inds = c(end.inds, toInd(x = df + dd, y = df, z = s, 
                                       x.max = n, y.max = n))
        }
      }
    }
    
    if(verbose) {
      message('Pre-computing bridging weights')
    }
    
    pBridge.pre = vector('list', N.max)
    
    for(i in 1:N.max){
      if(i==1) {
        # construct degenerate prob. of ending at a target node in 0 transitions
        pBridge.pre[[1]] = sparseMatrix(i = end.inds, j = end.inds, 
                                        x = rep(1, length(end.inds)), 
                                        dims = rep(n^2 * 3, 2))
      } else {
        pBridge.pre[[i]] = tx.mat %*% pBridge.pre[[i-1]]
      }
    }
    
    # initialize bridged path output
    paths.out = vector('list', M)
    
    # sample paths
    for(i in 1:M) {
      
      # set starting index for trajectory
      init.ind = toInd(x = ifelse(is.null(d0.last[i]), n, d0.last[i]), 
                       y = d0, z = s0[i], x.max = n, y.max = n)
      
      # extract path length
      N.i = N[i]
      
      if(verbose) {
        message(paste('Sampling', N.i, 'potential transitions for trajectory', 
                      i, sep = ' '))
      }
      
      # initialize trajectory and sample path
      current.ind = init.ind
      path.inds = numeric(N.i)
      
      # sample path transitions
      if(N.i == 0) {
        
        # extract raw path, including initial state
        path.inds = init.ind
        
        t.thick = c(t0, tf)
        
      } else {
        for(k in 1:N.i) {
          
          # define transition support
          out.inds = out.inds.lookup[[current.ind]]
          
          # extract bridging weights
          pBr = pBridge.pre[[N.i - k + 1]]
          if(length(out.inds) > 1) {
            wts = rowSums(pBr[out.inds,end.inds])
          } else {
            wts = sum(pBr[out.inds,end.inds])
          }
          
          # compute transition probability
          tx.dens = tx.mat[current.ind, out.inds] * wts
          tx.dens = tx.dens / sum(tx.dens)
          
          #
          # sample and update path
          #
          
          if(cond.sim & i==1) {
            current.ind = toInd(x = trajectory.conditional$depths[k],  
                                y = trajectory.conditional$depths[k+1], 
                                z = trajectory.conditional$stages[k+1], 
                                x.max = n, y.max = n)
          } else {
            if(length(out.inds) > 1) {
              current.ind = sample(x = out.inds, size = 1, prob = tx.dens)
            } else {
              current.ind = out.inds
            }
          }
          
          # update log-density
          ld[i] = ld[i] + log(tx.dens[which(current.ind == out.inds)])
          
          path.inds[k] = current.ind
          
        }
        
        if(cond.sim & i==1) {
          t.thick = trajectory.conditional$times
        } else {
          # sample arrival times, including initial time
          t.thick = c(t0, t0 + T.win * sort(runif(n = N[i])))
        }
      }
      
      #
      # package path
      #

      # extract raw path, including initial state
      path.full = cbind( 
        c(ifelse(is.null(d0.last[i]), n, d0.last[i]), d0, s0[i]), 
        sapply(path.inds, function(ind){
          fromInd(ind = ind, x.max = n, y.max = n)
        })
      )
      
      # remove self-transitions, and include initial state
      true.tx.inds = c(TRUE, diff(path.full[2,]) != 0)
      path.full = path.full[, true.tx.inds, drop = FALSE]
      
      # package trajectory
      p = list(
        depths = path.full[2,],
        stages = path.full[3,],
        times = t.thick[true.tx.inds],
        durations = diff(t.thick[true.tx.inds]),
        ld = ld[i]
      )
      
      if(ncol(path.full)==1) {
        p$durations = tf - t0
      }
      
      # compute density under true model 
      durations.tmp = c(p$durations, tf - p$times[length(p$times)])
      p$ld.true = dsdive.ld(depths = p$depths, durations = durations.tmp, 
                            times = p$times, stages = p$stages, beta = beta, 
                            lambda = lambda, sub.tx = sub.tx, surf.tx = surf.tx, 
                            depth.bins = depth.bins, t0.dive = t0.dive, 
                            d0.last = d0.last[i])
      
      # save trajectory
      paths.out[i] = list(p)
      
      class(paths.out[i]) = 'dsdive'
      
    }
  
  # package and return results
  paths.out
}