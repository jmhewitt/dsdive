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
#' @example examples/dsdive.impute_segments_cond.R
#'
#' @importFrom stats runif
#' @importFrom extraDistr rtpois dtpois
#' @import Matrix
#'
#' @export
#'
#'
dsdive.impute_segments_cond = function(depth.bins, depths, times, beta, 
                                       lambda, inflation.factor.lambda = 1.1,
                                       verbose = FALSE, t.sbreaks, t.cond) {

  if(!is.null(t.sbreaks)) {
    if(any(t.sbreaks <= times[1])) {
      stop('stage tx times cannot occur at or before first imputation time')
    }
    if(any(t.sbreaks > times[length(times)])) {
      stop('stage tx times cannot occur after last observation')
    }
  }
  
  n.bins = nrow(depth.bins)
  
  # we are imputing a whole dive
  s.range = 1:3
  
  # determine rate for uniformized poisson process
  rate.unif = inflation.factor.lambda * 
    max(outer(lambda[s.range], 2 * depth.bins[,2], '/'))
  
  # compute uniformized depth bin transition matrices
  tx.mat = lapply(s.range, function(s) {
    dsdive.tx.matrix.uniformized(depth.bins = depth.bins, beta = beta, 
                                 lambda = lambda, s0 = s, 
                                 rate.uniformized = rate.unif)
  })
  
  # transition matrix to force stage transitions
  tx.between = stage_tx.matrix.uniformized(m.list = tx.mat)
  
  # transition matrix to only allow transitions within stages
  tx.within = bdiag(tx.mat)
  
  # get dive window
  T.win = times[length(times)] - times[1]
  
  # sample number of virtual arrivals, their arrival times, and stages
  N = rpois(n = 1, lambda = rate.unif * T.win)
  t.prop = times[1] + c(0, T.win * runif(N))
  
  # munge and label times
  t.all = sort(unique(c(t.prop, t.cond, t.sbreaks)))
  t.sbreaks.ind = which(t.all %in% t.sbreaks)
  
  # get stages for all times
  stages.all = s.range[findInterval(t.all, t.sbreaks) + 1]
  
  #
  # build likelihood and transition matrices for all timepoints
  #
  
  k = nrow(tx.between)
  
  # default entry is no information about the depth bin
  L = matrix(1/k, nrow = k, ncol = length(t.all))
  
  # default transition matrix is within-stage
  B = lapply(t.all, function(s) tx.within)
  
  # map observations to likelihood matrix
  depths.windows = findInterval(times, t.all)
  for(i in 1:length(depths)) {
    # which virtual timepoint is this depth observation associated with?
    virt.ind = depths.windows[i]
    # what is the stage-aligned index in the transition matrix for this depth?
    obs.ind = toInd(x = 1, y = depths[i], x.max = 1, y.max = n.bins,
                    z = which(s.range == stages.all[virt.ind]))
    # specify that the depth bin is known at the observation point
    L[obs.ind, virt.ind] = 1
    L[-obs.ind, virt.ind] = 0
  }
  
  # map transition information
  for(virt.ind in t.sbreaks.ind) {
    if(virt.ind == 1) {
      stop("Invalid stage break index")
    }
    # encode that a new stage has entered by time virt.ind
    last.old.ind = n.bins * stages.all[virt.ind-1]
    new.inds = (last.old.ind+1):nrow(tx.between)
    L[1:last.old.ind, virt.ind] = 0
    L[new.inds, virt.ind] = 1/length(new.inds)
    # transition matrix to force stage break
    B[[virt.ind-1]] = tx.between
  }
  

  # sample path
  d.prop = ffbs.segment(B = B, L = L)

  # decode stage-break indices
  d.prop = sapply(d.prop, function(ind) {
    fromInd(ind = ind, x.max = 1, y.max = n.bins)[2]
  })
  
  
  #
  # package path segment
  #
  
  # assemble segment imputations
  path.full = d.prop
  times.full = t.all
  stages.full = stages.all
  
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
