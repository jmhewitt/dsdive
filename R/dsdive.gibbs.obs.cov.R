#' Gibbs sampler for parameters of a model for dives across discrete depth bins
#'
#' This function differs from \code{dsdive.gibbs.obs} in that this function
#' allows the model to be estimated with dive-specific covariates.
#'
#' @param dsobs.list list of \code{dsobs} objects, which describe the
#'   observation times and depths of a collection of dives
#' @param covs matrix of covariates associated with \code{dsobs.list}.  Each
#'   row ov \code{covs} should contain all covariates for a single dive.
#' @param t.stages.list list of initial stage transition times for dives
#'   observed in \code{dsobs.list}
#' @param beta.init Initial values for directional preference model parameters.
#'   See \code{dsdive.tx.params} for more details.  Should be a list, in which
#'   each component contains the initial value for coefficients for
#'   logit(pi^{(s)}.)
#' @param alpha.init Initial values for diving rate model parameters.  See
#'   \code{dsdive.tx.params} for more details.  Should be a list, in which
#'   each component contains the initial value for coefficients for
#'   log(lambda^{(s)}).
#' @param verbose \code{TRUE} to output sampler status while running
#' @param maxit number of Gibbs iterations to run
#' @param checkpoint.fn User-defined function to run during a checkpoint step;
#'   gives the user an opportunity to save partial output from the sampler
#' @param checkpoint.interval Number of seconds between calls to
#'   \code{checkpoint.fn}
#' @param beta1.prior List containing \code{mean} and \code{sd} parameters for
#'   independent Normal prior distributions on the coefficients for the
#'   dive-stage model parameter \eqn{\pi^{(1)}}.  The length of the \code{mean}
#'   and \code{sd} vectors should be equal to the number of coefficients for
#'   this model component.
#' @param beta2.prior List containing \code{mean} and \code{sd} parameters for
#'   independent Normal prior distributions on the coefficients for the
#'   dive-stage model parameter \eqn{\pi^{(3)}}.  The length of the \code{mean}
#'   and \code{sd} vectors should be equal to the number of coefficients for
#'   this model component.
#' @param alpha1.prior List containing \code{mean} and \code{sd} parameters for
#'   independent Normal prior distributions on the coefficients for the
#'   dive-stage model parameter \eqn{\lambda^{(1)}}.  The length of the
#'   \code{mean} and \code{sd} vectors should be equal to the number of
#'   coefficients for this model component.
#' @param alpha2.prior List containing \code{mean} and \code{sd} parameters for
#'   independent Normal prior distributions on the coefficients for the
#'   dive-stage model parameter \eqn{\lambda^{(2)}}.  The length of the
#'   \code{mean} and \code{sd} vectors should be equal to the number of
#'   coefficients for this model component.
#' @param alpha3.prior List containing \code{mean} and \code{sd} parameters for
#'   independent Normal prior distributions on the coefficients for the
#'   dive-stage model parameter \eqn{\lambda^{(3)}}.  The length of the
#'   \code{mean} and \code{sd} vectors should be equal to the number of
#'   coefficients for this model component.
#' @param tstep Time between observations in \code{dsobs.list}
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first
#'   column defines the depth at the center of each depth bin, and the second
#'   column defines the half-width of each bin.
#' @param T1.prior.params \code{shape} and \code{rate} parameters for Gamma
#'   prior on the descent-stage duration.
#' @param T2.prior.params \code{shape} and \code{rate} parameters for Gamma
#'   prior on the bottom-stage duration.
#' @param max.width The stage transition times are updated with a piecewise
#'   proposal distribution.  \code{max.width} controls the maximum width of the
#'   intervals for the proposal distribution.  This is a tuning parameter that
#'   controls the numerical stability of the proposal distribution, which is
#'   sampled via inverse CDF techniques.
#' @param max.width.offset The t0 and tf offsets are updated with a piecewise
#'   proposal distribution.  \code{max.width.offset} controls the maximum width
#'   of the intervals for the proposal distribution.  This is a tuning parameter
#'   that controls the numerical stability of the proposal distribution, which
#'   is sampled via inverse CDF techniques.
#' @param t0.prior.params \code{shape1} and \code{shape2} parameters for the
#'   scaled and shifted Beta prior distribution for the t0 offset.
#' @param tf.prior.params \code{shape1} and \code{shape2} parameters for the
#'   scaled and shifted Beta prior distribution for the tf offset.
#' @param offsets vector with initial values for t0 offsets.
#' @param offsets.tf vector with initial values for tf offsets.
#' @param pi.formula List of formula objects, each of which defines the linear
#'   model that combines covariates in the \code{covs} matrix into
#'   dive-specific diving preferences logit(pi^{(s)}).
#' @param lambda.formula List of formula objects, each of which defines the
#'   linear model that combines covariates in the \code{covs} matrix into
#'   dive-specific diving rates log(lambda^{(s)}).
#' @param cl Shared-memory cluster to be used to distribute some computations.
#'   The cluster must be fully initialized before running this function.
#'   The cluster requires random seeds and \code{Rdsm}-initialization.
#' @param warmup number of iterations during which the proposal distributions 
#'   will be updated at each step
#' @param delta If \code{delta>0}, then the observation matrix and raw 
#'   components computed will be for a transition matrix whose generator is 
#'   perturbed to allow much faster computation.
#'   
#' @example examples/dsdive.gibbs.obs.cov.R
#'
#' @importFrom parallel clusterEvalQ
#' @importFrom stats model.matrix
#' @import Rdsm
#'
#' @export
#'
dsdive.gibbs.obs.cov = function(
  dsobs.list, covs, t.stages.list, beta.init, alpha.init, verbose = FALSE,
  maxit, checkpoint.fn, checkpoint.interval = 3600, beta1.prior, beta2.prior,
  alpha1.prior, alpha2.prior, alpha3.prior, tstep, depth.bins,
  T1.prior.params, T2.prior.params, max.width, max.width.offset,
  t0.prior.params, tf.prior.params, offsets, offsets.tf, cl,
  pi.formula, lambda.formula, warmup = Inf, delta = 1e-10) {


  n = length(dsobs.list)

  alpha.priors.list = list(alpha1.prior, alpha2.prior, alpha3.prior)
  beta.priors.list = list(beta1.prior, beta2.prior)


  #
  # expand covariate design matrices
  #

  if(!inherits(pi.formula, 'list')) {
    pi.formula = list(pi.formula, pi.formula)
  }

  if(!inherits(lambda.formula, 'list')) {
    lambda.formula = list(lambda.formula, lambda.formula, lambda.formula)
  }

  pi.designs = lapply(pi.formula, function(f) model.matrix(f, covs))
  lambda.designs = lapply(lambda.formula, function(f) model.matrix(f, covs))


  #
  # initialize nodes and shared memory pass-throughs
  #

  shared.env = gibbs_init_shared(cl = cl, envir = environment())


  #
  # initialize sampler output
  #

  theta = list(beta1 = shared.env$beta1[], beta2 = shared.env$beta2[],
               alpha1 = shared.env$alpha1[], alpha2 = shared.env$alpha2[],
               alpha3 = shared.env$alpha3[])

  trace = list(
    beta1 = matrix(NA, nrow = maxit, ncol = length(beta.init[[1]])),
    beta2 = matrix(NA, nrow = maxit, ncol = length(beta.init[[2]])),
    alpha1 = matrix(NA, nrow = maxit, ncol = length(alpha.init[[1]])),
    alpha2 = matrix(NA, nrow = maxit, ncol = length(alpha.init[[2]])),
    alpha3 = matrix(NA, nrow = maxit, ncol = length(alpha.init[[3]]))
  )
  trace.t.stages = vector('list', length = maxit)
  trace.offsets = matrix(0, nrow = maxit, ncol = n)
  trace.offsets.tf = matrix(0, nrow = maxit, ncol = n)

  tick.checkpoint = proc.time()[3]

  # proposal distribution caches
  proposaldists.theta = vector('list', 3)
  proposaldists.offsets = vector('list', n)
  proposaldists.offsets.tf = vector('list', n)


  #
  # gibbs sample
  #

  for(it in 1:maxit) {

    if(verbose) {
      message(paste('Iteration:', it, sep = ' '))
      tick = proc.time()[3]
    }


    #
    # sample model parameters from full conditional posterior densities
    #

    # update stage 1 parameters
    theta.raw = dsdive.obs.sampleparams_shared(
      s0 = 1, theta = theta, alpha.priors.list = alpha.priors.list,
      beta.priors.list = beta.priors.list, cl = cl, shared.env = shared.env)

    theta = theta.raw$theta

    # update stage 1 tx matrices on nodes
    if(theta.raw$accepted) {
      clusterEvalQ(cl, {
        s = 1
        for(local_ind in 1:local_n) {
          id = local_ids[local_ind]
          P.raw[[local_ind]][[s]] = dsdive.obstxmat.cov(
            pi.designs = pi.designs, lambda.designs = lambda.designs,
            beta1 = beta1[], beta2 = beta2[], alpha1 = alpha1[],
            alpha2 = alpha2[], alpha3 = alpha3[], s0 = s,
            ind = id, tstep = tstep, include.raw = TRUE,
            depth.bins = depth.bins, delta = delta)
        }
      })
    }


    # update stage 2 parameters
    theta.raw = dsdive.obs.sampleparams_shared(
      s0 = 2, theta = theta, alpha.priors.list = alpha.priors.list,
      beta.priors.list = beta.priors.list, cl = cl, shared.env = shared.env)

    theta = theta.raw$theta

    # update stage 2 tx matrices on nodes
    if(theta.raw$accepted) {
      clusterEvalQ(cl, {
        s = 2
        for(local_ind in 1:local_n) {
          id = local_ids[local_ind]
          P.raw[[local_ind]][[s]] = dsdive.obstxmat.cov(
            pi.designs = pi.designs, lambda.designs = lambda.designs,
            beta1 = beta1[], beta2 = beta2[], alpha1 = alpha1[],
            alpha2 = alpha2[], alpha3 = alpha3[], s0 = s,
            ind = id, tstep = tstep, include.raw = TRUE,
            depth.bins = depth.bins, delta = delta)
        }
      })
    }


    # update stage 3 parameters
    theta.raw = dsdive.obs.sampleparams_shared(
      s0 = 3, theta = theta, alpha.priors.list = alpha.priors.list,
      beta.priors.list = beta.priors.list, cl = cl, shared.env = shared.env)

    theta = theta.raw$theta

    # update stage 3 tx matrices on nodes
    if(theta.raw$accepted) {
      clusterEvalQ(cl, {
        s = 3
        for(local_ind in 1:local_n) {
          id = local_ids[local_ind]
          P.raw[[local_ind]][[s]] = dsdive.obstxmat.cov(
            pi.designs = pi.designs, lambda.designs = lambda.designs,
            beta1 = beta1[], beta2 = beta2[], alpha1 = alpha1[],
            alpha2 = alpha2[], alpha3 = alpha3[], s0 = s,
            ind = id, tstep = tstep, include.raw = TRUE,
            depth.bins = depth.bins, delta = delta)
        }
      })
    }

    if(verbose) {
      print(theta)
    }


    #
    # update stage transition time parameters and dive offsets
    #

    clusterEvalQ(cl, {
      for(local_ind in 1:local_n) {

        id = local_ids[local_ind]

        # sample dive stage transition times
        d = dsobs.aligned[[id]]
        t.stages[id,] = dsdive.obs.sample.stages(
          depths = d$depths, times = d$times, t.stages = t.stages[id,],
          P.raw = P.raw[[local_ind]],
          T.range = c(d$times[1], d$times[length(d$times)]),
          depth.bins = depth.bins, T1.prior.params = T1.prior.params,
          T2.prior.params = T2.prior.params, max.width = max.width,
          debug = FALSE)$t.stages

        # sample dive start offset
        if(!is.null(t0.prior.params)) {
          d0 = dsdive.obs.sample.offsets(
            dsobs.aligned = d, dsobs.unaligned = dsobs.list[[id]],
            offset = offsets[id], offset.tf = offsets.tf[id],
            t.stages = t.stages[id,], P.raw = P.raw[[local_ind]],
            depth.bins = depth.bins, tstep = tstep,
            max.width = max.width.offset, prior.params = t0.prior.params,
            sample.start = TRUE, lpapprox = NULL, output.lpapprox = TRUE)
          # extract new offset
          dsobs.aligned[[id]] = d0$dsobs.aligned
          offsets[id] = d0$offset
        }

        # sample dive end offset
        if(!is.null(tf.prior.params)) {
          d0 = dsdive.obs.sample.offsets(
            dsobs.aligned = d, dsobs.unaligned = dsobs.list[[id]],
            offset = offsets[id], offset.tf = offsets.tf[id],
            t.stages = t.stages[id,], P.raw = P.raw[[local_ind]],
            depth.bins = depth.bins, tstep = tstep,
            max.width = max.width.offset, prior.params = tf.prior.params,
            sample.start = FALSE, lpapprox = NULL, output.lpapprox = TRUE)
          # extract new offset
          dsobs.aligned[[id]] = d0$dsobs.aligned
          offsets.tf[id] = d0$offset
        }

      }
    })


    #
    # save trace
    #

    trace$beta1[it,] = theta$beta1
    trace$beta2[it,] = theta$beta2
    trace$alpha1[it,] = theta$alpha1
    trace$alpha2[it,] = theta$alpha2
    trace$alpha3[it,] = theta$alpha3
    trace.t.stages[[it]] = shared.env$t.stages[]
    trace.offsets[it,] = shared.env$offsets[]
    trace.offsets.tf[it,] = shared.env$offsets.tf[]

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
      checkpoint.fn(list(theta = trace[1:it],
                         trace.t.stages = trace.t.stages[1:it],
                         trace.offsets = trace.offsets[1:it,],
                         trace.offsets.tf = trace.offsets.tf[1:it,]))
      tick.checkpoint = proc.time()[3]
    }

  }

  list(
    theta = trace,
    trace.t.stages = trace.t.stages,
    trace.offsets = trace.offsets,
    trace.offsets.tf = trace.offsets.tf
  )
}
