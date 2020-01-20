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
    if(any(t.sbreaks < times[1])) {
      stop('t.sbreaks cannot occur before first imputation time')
    }
  }
  
  # end stage for segment
  sf = s0 + ifelse(is.null(t.sbreaks), 0 , 
                   sum(t.sbreaks < times[length(times)]))
  
  # sampling requirements for each segment
  T.win = diff(times)
  min.tx = abs(diff(depths))
  
  # determine rate for uniformized poisson process
  rate.unif = inflation.factor.lambda * 
    max(outer(lambda[s0:sf], 2 * depth.bins[,2], '/'))
  
  # compute uniformized stage transition matrices
  tx.mat = lapply(s0:sf, function(s) {
    dsdive.tx.matrix.uniformized(depth.bins = depth.bins, beta = beta, 
                                 lambda = lambda, s0 = s, 
                                 rate.uniformized = rate.unif)
  })
  
  # sample dive segments
  t.prop = vector('list', length(T.win))
  d.prop = vector('list', length(T.win))
  stages.prop = vector('list', length(T.win))
  for(i in 1:length(T.win)) {
    
    # sample number of pseudo-arrivals
    if(method.N == 'truncpois') {
      N = rtpois(n = 1, lambda = rate.unif * T.win[i], a = min.tx[i])
    } else if(method.N == 'exact') {
      dN = dN.bridged(B = tx.mat, x0 = depths[i], xN = depths[i+1], 
                      N.max = N.max, rate.uniformized = rate.unif, 
                      t = T.win[i], log = TRUE)
      N = sample.gumbeltrick(dN) - 1
    }
    
    # if it is in the interval, add t.sbreak as a transition time as well, since this is true!
    
    # sample arrival times and stages
    t.prop[[i]] = times[i] + c(0, T.win[i] * sort(runif(N)))
    stages.prop[[i]] = (s0:sf)[findInterval(t.prop[[i]], t.sbreaks) + 1]
    
    # pair transition matrices with each timepoint
    B = lapply(stages.prop[[i]], function(s) tx.mat[[s-s0+1]])
    
    
    # sample path
    d.prop[[i]] = ffbs.segment(B = B, x0 = depths[i], xN = depths[i+1], N = N)
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
