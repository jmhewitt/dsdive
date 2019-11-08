#' Impute dive trajectories that meet start and end times
#'
#' Sample a dive trajectory such that its start and end locations and times are 
#' fixed.  Sampling uses a thinned Poisson process and bridged Markov chain 
#' transition probabilities.  Sampling is not from the model's theoretical 
#' imputation distribution, but an approximation thereof.  The trajectory's
#' log-density is returned with the output.
#' 
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
#' @param precompute.bridges If \code{TRUE}, then the bridged transition 
#'   matrices will be precomputed.  Enabling this option will increase the 
#'   memory overhead of the method, but will reduce its runtime.
#' @param lambda.max Arrival rate for the parent Poisson process that will
#'   be thinned.  \code{lambda.max} will be scaled by 
#'   \code{inflation.factor.lambda}, and if \code{lambda.max==NULL} then the 
#'   method will compute this on its own.
#' 
#' @example examples/dsdive.bridgesample.R
#' 
#' @importFrom extraDistr rtpois dtpois
#' @import Matrix
#' 
#' @export
#' 
#'
dsdive.bridgesample = function(depth.bins, d0, d0.last, df, beta, lambda, 
                               sub.tx, surf.tx, t0, tf, s0, 
                               inflation.factor.lambda = 1.1, verbose = FALSE,
                               precompute.bridges = TRUE, lambda.max = NULL) {
  
  #
  # build basic simulation parameters
  #
  
  # initialize log-density for sample
  ld = 0
  
  # transition time window
  T.win = tf - t0
  
  # minimum number of transitions required for bridging
  min.tx = abs(df-d0)
  
  
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
  N = rtpois(n = 1, lambda = lambda.tmp, a = min.tx)
  
  if(verbose) {
    message(paste('Simulating', N, 'potential transitions', sep =' '))
  }
  
  # update log-density for sample
  ld = ld + dtpois(x = N, lambda = lambda.tmp, a = min.tx, log = TRUE) + 
    lfactorial(N) - N * log(T.win)
  
  # sample arrival times, including initial time
  t.thick = c(t0, t0 + T.win * sort(runif(n = N)))
  
  
  #
  # bridge sample a trajectory
  #
  
  if(verbose) {
    message('Computing transition matrices at potential transition times')
  }
  
  # build complete transition matrices for each arrival time
  min.depth = max(1, d0 - N)
  max.depth = min(nrow(depth.bins), d0 + N)
  
  # compute transition matrices at times t(0), t(1), ... t(N-1)
  tx.mat = vector('list', N)
  for(i in 1:N) {
    tx.mat[[i]] = dsdive.tx.matrix(t0 = t.thick[i], depth.bins = depth.bins, 
                                   beta = beta, lambda = lambda, 
                                   sub.tx = sub.tx, surf.tx = surf.tx,
                                   inflation.factor.lambda = 
                                     inflation.factor.lambda,
                                   min.depth = min.depth, max.depth = max.depth)
  }
  
  # add a "null" depth bin to allow trajectory initialization
  n = nrow(depth.bins) + 1
  
  # set starting index for trajectory
  current.ind = toInd(x = ifelse(is.null(d0.last), n, d0.last), y = d0, z = s0, 
                      x.max = n, y.max = n)
  
  # get ending indices for trajectory (any previous depth, any dive stage)
  end.inds = c()
  for(s in s0:3) {
    for(dd in c(-1,0,1)) {
      if(((df + dd) > 0) & ((df + dd) < n)) {
        end.inds = c(end.inds, toInd(x = df + dd, y = df, z = s, 
                                     x.max = n, y.max = n))
      }
    }
  }
  
  if(precompute.bridges) {
    
    if(verbose) {
      message('Pre-computing bridging weights')
    }
    
    pBridge.pre = vector('list', N)
  
    # construct degenerate prob. of ending at a target node in 0 transitions
    pBridge.pre[[N]] = sparseMatrix(i = end.inds, j = end.inds, 
                                    x = rep(1, length(end.inds)), 
                                    dims = rep(n^2 * 3, 2))
    
    if(N > 1) {
      pBridge.pre[[N-1]] = tx.mat[[N]]
    }
    
    if(N-2 > 0) {
      for(k in (N-2):1) {
        pBridge.pre[[k]] = tx.mat[[k+1]] %*% pBridge.pre[[k+1]]
      }
    }
  }
  
  # initialize bridged path
  path.inds = numeric(length = N)
  
  
  # sample transitions
  for(k in 1:N) {
    
    if(verbose) {
      message(paste('Sampling', k, 'of', N, 'potential transitions', sep = ' '))
    }
    
    #
    # compute local transition probabilities
    #
    
    if(precompute.bridges) {
      pBridge = pBridge.pre[[k]]
    } else {
      # compute bridging probabilities 
      pBridge = tx.mat[[k+1]]
      if(k + 2 <= N) {
        for(j in (k+2):N) {
          pBridge = pBridge %*% tx.mat[[j]]
        }
      } 
      # at last transition, bridging probabilities degenerate to indicator fns.
      else if(k==N) {
        pBridge = sparseMatrix(i = end.inds, j = end.inds, 
                               x = rep(1, length(end.inds)))
      }
    }
    
    # define transition support
    out.inds = which(tx.mat[[k]][current.ind,] > 0)
    
    # extract bridging weights
    if(length(out.inds) > 1) {
      wts = rowSums(pBridge[out.inds,end.inds])
    } else {
      wts = sum(pBridge[out.inds,end.inds])
    }
    
    # compute transition probability
    tx.dens = tx.mat[[k]][current.ind, out.inds] * wts
    tx.dens = tx.dens / sum(tx.dens)
    
    
    #
    # sample and update path
    #
    
    if(length(out.inds) > 1) {
      current.ind = sample(x = out.inds, size = 1, prob = tx.dens)
    } else {
      current.ind = out.inds
    }
    path.inds[k] = current.ind
    
    # update log-density
    ld = ld + log(tx.dens[which(current.ind == out.inds)])
    
  }

    
  #
  # package and return results
  #
  
  # extract raw path, including initial state
  path.full = cbind( 
    c(ifelse(is.null(d0.last), n, d0.last), d0, s0), 
    sapply(path.inds, function(ind){
      fromInd(ind = ind, x.max = n, y.max = n)
    })
  )
  
  # remove self-transitions, and include initial state
  true.tx.inds = c(TRUE, diff(path.full[2,]) != 0)
  path.full = path.full[, true.tx.inds, drop = FALSE]
  
  # package and return results
  res = list(
    depths = path.full[2,],
    stages = path.full[3,],
    t = t.thick[true.tx.inds],
    durations = diff(t.thick[true.tx.inds]),
    ld = ld
  )
  
  res
}