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
#' @example examples/dsdive.impute_segments.R
#'
#' @importFrom stats runif
#' @importFrom extraDistr rtpois dtpois
#' @import Matrix
#'
#' @export
#'
#'
dsdive.impute_segments = function(depth.bins, depths, times, beta, 
                                  lambda, s0, inflation.factor.lambda = 1.1,
                                  verbose = FALSE, method.N = 'exact',
                                  N.max = NULL, t.sbreaks = NULL) {

  if(!is.null(t.sbreaks)) {
    if(any(t.sbreaks <= times[1])) {
      stop('stage tx times cannot occur at or before first imputation time')
    }
    if(any(t.sbreaks > times[length(times)])) {
      stop('stage tx times cannot occur after last observation')
    }
  }
  
  # end stage for segment
  sf = s0 + length(t.sbreaks)
  s.range = s0:sf
  
  # sampling requirements for each segment
  T.win = diff(times)
  min.tx = abs(diff(depths))
  
  # determine rate for uniformized poisson process
  rate.unif = inflation.factor.lambda * 
    max(outer(lambda[s.range], 2 * depth.bins[,2], '/'))
  
  # compute uniformized depth bin transition matrices
  tx.mat = lapply(s.range, function(s) {
    dsdive.tx.matrix.uniformized(depth.bins = depth.bins, beta = beta, 
                                 lambda = lambda, s0 = s, 
                                 rate.uniformized = rate.unif)
  })
  
  # build within and between-stage transition matrices
  if(s0!=sf) {
    # transition matrix to force stage transitions
    tx.between = stage_tx.matrix.uniformized(m.list = tx.mat)
    # only allow transitions within stages
    tx.within = bdiag(tx.mat)
  }
  
  # extract number of depth bins
  n.bins = nrow(tx.mat[[1]])
  
  # sample dive segments
  t.prop = vector('list', length(T.win))
  d.prop = vector('list', length(T.win))
  stages.prop = vector('list', length(T.win))
  for(i in 1:length(T.win)) {
    
    # check to see if any stage breaks occur within interval
    breaks.now = which((times[i] < t.sbreaks) & (t.sbreaks <= times[i+1]))
    n.breaks = length(breaks.now)
    
    # adjust minimum number of transitions to accomodate stage txs.
    if(n.breaks > 0) {
      min.tx[i] = min.tx(d0 = depths[i], df = depths[i+1], ns = n.breaks)
    }
    
    # sample number of pseudo-arrivals
    if(method.N == 'truncpois') {
      N = rtpois(n = 1, lambda = .95 * rate.unif * T.win[i], a = min.tx[i])
    } else if(method.N == 'exact') {
      dN = dN.bridged(B = tx.mat[[1]], x0 = depths[i], xN = depths[i+1], 
                      N.max = N.max, rate.uniformized = rate.unif, 
                      t = T.win[i], log = TRUE)
      N = sample.gumbeltrick(dN) - 1
    }
    
    # sample arrival times
    t.prop[[i]] = unique(sort(c(
      t.sbreaks[breaks.now], 
      times[i] + c(0, T.win[i] * runif(N-n.breaks))
    )))
    
    # determine stages for each timepoint
    stages.prop[[i]] = s.range[findInterval(t.prop[[i]], t.sbreaks) + 1]
    
    # we end up with N+1 arrival times
    
    # reshuffle arrival times if necessary
    if(n.breaks > 0) {
      if(depths[i+1] == 1) {
        # need at least two stage 3 transitions to guarantee surfacing
        if(sum(stages.prop[[i]] == 3) < 2) {
          # redistribute number of arrivals between stage 3 transition and end
          t.s3 = t.prop[[i]][stages.prop[[i]] == 3]
          win.3 = times[i+1] - t.s3
          N.3 = rtpois(n = 1, lambda = .95 * rate.unif * win.3, a = 1, b = N-1)
          t.prop[[i]] = unique(sort(c(
            times[i] + c(0, (t.s3 - times[i]) * runif(N-N.3-1)),
            t.s3 + c(0, win.3 * runif(N.3))
          )))
          # reset the stages
          stages.prop[[i]] = s.range[findInterval(t.prop[[i]], t.sbreaks) + 1]
        }
      }
    }
    
    # build transition and likelihood matrices for each timepoint
    if(n.breaks > 0) {
      
      B = lapply(stages.prop[[i]], function(s) tx.within)
      stx.inds = which(diff(stages.prop[[i]])==1)
      for(j in stx.inds) {
        B[[j]] = tx.between
      }
      
      ind.start = toInd(x = 1, y = depths[i], x.max = 1, y.max = n.bins,
                        z = which(s.range == stages.prop[[i]][1]))
      ind.end = toInd(x = 1, y = depths[i+1], x.max = 1, y.max = n.bins, 
                      z = which(s.range == 
                                  stages.prop[[i]][length(stages.prop[[i]])]))
      
      L = matrix(0, nrow = nrow(tx.between), ncol = N+1)
      L[ind.start,1] = 1                  # fix start bin
      L[ind.end,N+1] = 1                  # fix end bin
      L[,-c(1,N+1)] = 1/nrow(tx.between)  # free transitions for other steps
      
    } else {
      B = lapply(stages.prop[[i]], function(s) tx.mat[[s-s0+1]])
      L = matrix(0, nrow = n.bins, ncol = N+1)
      L[depths[i],1] = 1        # fix start bin
      L[depths[i+1],N+1] = 1    # fix end bin
      L[,-c(1,N+1)] = 1/n.bins  # free transitions for intermediate bins
    }
      
    # sample path
    d.prop[[i]] = ffbs.segment(B = B, L = L)
    
    # decode stage-break indices
    if(n.breaks > 0) {
      d.prop[[i]] = sapply(d.prop[[i]], function(ind) 
        fromInd(ind = ind, x.max = 1, y.max = n.bins)[2]
      )
    }
  }
  
  
  #
  # package path segment
  #
  
  # assemble segment imputations
  path.full = do.call(c, d.prop)
  times.full = do.call(c, t.prop) 
  stages.full = do.call(c, stages.prop) 
  
  # remove self-transitions, and include initial state
  true.tx.inds = c(TRUE, diff(path.full) != 0)
  path.full = path.full[true.tx.inds]
  times.full = times.full[true.tx.inds]
  stages.full = stages.full[true.tx.inds]
  
  res = list(
    depths = path.full,
    durations = diff(times.full),
    times = times.full,
    stages = stages.full
  )
  
  class(res) = 'dsdive'
  
  res
}
