#' Impute a sequence of dive segments from a fixed transition matrix
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
#' @param t.stage2 vector of times at which stage 2 was entered.
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the
#'   method used to determine stage transition probability curves
#' @param ld.compute \code{TRUE} to compute likelihood values as well.
#'
#' @example examples/dsdive.batchedbridge.R
#'
#' @importFrom stats runif
#' @importFrom extraDistr rtpois dtpois
#' @import Matrix
#'
#' @export
#'
#'
dsdive.batchedbridge = function(depth.bins, depths, d0.last, times, tx.max,
                                beta, lambda, sub.tx, surf.tx, s0, sf,
                                inflation.factor.lambda = 1.1,
                                verbose = FALSE,
                                lambda.max = NULL, t0.dive, t.stages,
                                model) {

  s.range = s0:sf
  t0 = times[1]

  # set stage 2 entry time if transition occured before t0
  if(t0 >= t.stages[1]) {
    t.stage2 = t.stages[1]
  } else {
    t.stage2 = NA
  }


  #
  # construct homogeneous poisson process for thinning
  #

  # get rate for thick poisson process
  if(is.null(lambda.max)) {
    lambda.max = max(outer(lambda, 2 * depth.bins[,2], '/'))
  }
  lambda.thick = inflation.factor.lambda * lambda.max


  # compute the base transition matrix
  #
  # note: in ideal cases, the transition matrix only has a weak dependence on
  # t0, so we accelerate the sampling process by basing all proposals off of a
  # single transition matrix
  tx.mat.raw = dsdive.tx.matrix.simplified(t0 = t0, depth.bins = depth.bins, 
                                           beta = beta, lambda = lambda, 
                                           sub.tx = sub.tx, surf.tx = surf.tx, 
                                           s0 = s0, sf = sf, 
                                           inflation.factor.lambda = 1, 
                                           t0.dive = t0.dive, 
                                           lambda.max = lambda.thick, 
                                           t.stage2 = t.stage2, model = model)
  tx.mat = tx.mat.raw$m
  out.inds.lookup = tx.mat.raw$out.inds

  # add a "null" depth bin to allow trajectory initialization
  n = nrow(depth.bins) + 1

  # support and density component for number of transitions
  n.range = 0:tx.max
  n.lfact = lfactorial(n.range)

  # compute n-ahead transition matrices
  Rn = vector('list', length(n.range))
  Rn[[1]] = .sparseDiagonal(n = nrow(tx.mat))
  for(i in 2:length(n.range)) {
    Rn[[i]] = Rn[[i-1]] %*% tx.mat
  }

  # initialize output
  path.out = list(
    depths = depths[1], times = times[1], stages = s0, durations = NULL
  )

  # bridge segments
  nt = length(depths) - 1
  for(i in 1:nt) {

    # extract segment start/end points and time
    df = depths[i+1]
    t0 = times[i]
    tf = times[i+1]
    dt = tf - t0
    d0 = depths[i]
    
    if(df == 1) {
      df = n
    }

    # get ending indices for segment (any dive stage)
    end.inds = c()
    for(s in s.range) {
      end.inds = c(end.inds, toInd(x = 1, y = df,
                                   z = which(s==s.range),
                                   x.max = 1, y.max = n))
    }

    # set starting index for trajectory
    start.ind = toInd(x = 1, y = d0, z = 1, x.max = 1, y.max = n)

    # un-normalized (partial) pmf for number of transitions
    mT = lambda.thick * dt
    dN.log = - mT + n.range * log(mT) - n.lfact
    for(j in 1:length(n.range)) {
      dN.log[j] = dN.log[j] + log(sum(Rn[[j]][start.ind, end.inds]))
    }

    # sample number of transitions for segment
    N = sample.gumbeltrick(log.p = dN.log)

    # initialize trajectory and sample path
    current.ind = start.ind
    path.inds = numeric(N)

    if(N==0) {
      path.inds = init.ind
      t.thick = c(t0, tf)
    } else {
      for (k in 1:N) {
  
        # define transition support
        out.inds = out.inds.lookup[[current.ind]]
  
        # extract bridging weights
        if(length(out.inds) > 1) {
          wts = rowSums(Rn[[N-k+1]][out.inds, end.inds, drop = FALSE])
        } else {
          wts = sum(Rn[[N-k+1]][out.inds, end.inds])
        }
  
        #  compute transition probability
        tx.dens = tx.mat[current.ind, out.inds] * wts
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

      }
    }
    
    # sample arrival times, including initial time
    t.thick = c(t0, t0 + dt * sort(runif(n = N)))

    #
    # package path segment
    #

    # extract raw path, including initial state
      path.full = cbind(
        c(1, d0, 1),
        sapply(path.inds, function(ind){
          fromInd(ind = ind, x.max = 1, y.max = n)
        })
      )
      path.full[3,] = s.range[path.full[3,]]

    # remove self-transitions, and include initial state
    true.tx.inds = c(TRUE, diff(path.full[2,]) != 0)
    path.full = path.full[, true.tx.inds, drop = FALSE]
    
    # package path segment
    path.out$depths = c(path.out$depths, path.full[2,-1])
    path.out$stages = c(path.out$stages, path.full[3,-1])
    path.out$times = c(path.out$times, t.thick[true.tx.inds][-1])
      
  }
  
  # compute durations
  if(length(path.out$depths)==1) {
    path.out$durations = diff(range(times))
  } else {
    path.out$durations = diff(path.out$times)
  }
  
  class(path.out) = 'dsdive'
  
  path.out
}
