% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsdive.gibbs.obs.R
\name{dsdive.gibbs.obs}
\alias{dsdive.gibbs.obs}
\title{Gibbs sampler for parameters of a model for dives across discrete depth bins}
\usage{
dsdive.gibbs.obs(
  dsobs.list,
  t.stages.list,
  beta.init,
  lambda.init,
  verbose = FALSE,
  maxit,
  checkpoint.fn,
  checkpoint.interval = 3600,
  pi1.prior,
  pi2.prior,
  lambda1.prior,
  lambda2.prior,
  lambda3.prior,
  tstep,
  depth.bins,
  T1.prior.params,
  T2.prior.params,
  max.width,
  max.width.offset,
  t0.prior.params,
  tf.prior.params,
  offsets,
  offsets.tf,
  warmup = Inf,
  cl = NULL,
  delta
)
}
\arguments{
\item{dsobs.list}{list of \code{dsobs} objects, which describe the 
observation times and depths of a collection of dives}

\item{t.stages.list}{list of initial stage transition times for dives 
observed in \code{dsobs.list}}

\item{beta.init}{Initial values for directional preference model parameters.
See \code{dsdive.tx.params} for more details.}

\item{lambda.init}{Initial values for diving rate model parameters.  See 
\code{dsdive.tx.params} for more details.}

\item{verbose}{\code{TRUE} to output sampler status while running}

\item{maxit}{number of Gibbs iterations to run}

\item{checkpoint.fn}{User-defined function to run during a checkpoint step;
gives the user an opportunity to save partial output from the sampler}

\item{checkpoint.interval}{Number of seconds between calls to
\code{checkpoint.fn}}

\item{pi1.prior}{\code{shape1} and \code{shape2} parameters for Beta prior 
distribution on the dive-stage model parameter \eqn{\pi^{(1)}}}

\item{pi2.prior}{\code{shape1} and \code{shape2} parameters for Beta prior 
distribution on the ascent-stage model parameter \eqn{\pi^{(3)}}.  The 
notation is a little odd because this is the SECOND preference parameter the 
model estimates.}

\item{lambda1.prior}{\code{shape} and \code{rate} parameters for Gamma prior 
distribution on the descent-stage diving rate \eqn{\lambda^{(1)}}.}

\item{lambda2.prior}{\code{shape} and \code{rate} parameters for Gamma prior 
distribution on the bottom-stage diving rate \eqn{\lambda^{(2)}}.}

\item{lambda3.prior}{\code{shape} and \code{rate} parameters for Gamma prior 
distribution on the ascent-stage diving rate \eqn{\lambda^{(3)}}.}

\item{tstep}{Time between observations in \code{dsobs.list}}

\item{depth.bins}{\eqn{n x 2} Matrix that defines the depth bins.  The first 
column defines the depth at the center of each depth bin, and the second 
column defines the half-width of each bin.}

\item{T1.prior.params}{\code{shape} and \code{rate} parameters for Gamma 
prior on the descent-stage duration.}

\item{T2.prior.params}{\code{shape} and \code{rate} parameters for Gamma 
prior on the bottom-stage duration.}

\item{max.width}{The stage transition times are updated with a piecewise 
proposal distribution.  \code{max.width} controls the maximum width of the 
intervals for the proposal distribution.  This is a tuning parameter that 
controls the numerical stability of the proposal distribution, which is 
sampled via inverse CDF techniques.}

\item{max.width.offset}{The t0 and tf offsets are updated with a piecewise 
proposal distribution.  \code{max.width.offset} controls the maximum width 
of the intervals for the proposal distribution.  This is a tuning parameter 
that controls the numerical stability of the proposal distribution, which 
is sampled via inverse CDF techniques.}

\item{t0.prior.params}{\code{shape1} and \code{shape2} parameters for the 
scaled and shifted Beta prior distribution for the t0 offset.}

\item{tf.prior.params}{\code{shape1} and \code{shape2} parameters for the 
scaled and shifted Beta prior distribution for the tf offset.}

\item{offsets}{vector with initial values for t0 offsets.}

\item{offsets.tf}{vector with initial values for tf offsets.}

\item{warmup}{number of iterations during which the proposal distributions 
will be updated at each step}

\item{cl}{cluster to be used to distribute some computations}

\item{delta}{If \code{delta>0}, then the probability transition matrices
computed will use a transition matrix whose generator is 
perturbed to allow much faster computation.  See \code{dsdive.obstx.matrix}
for more details.}
}
\description{
Gibbs sampler for parameters of a model for dives across discrete depth bins
}
\examples{
data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

t.stages = sim$times[c(FALSE,diff(sim$stages)==1)]
  
tstep = diff(sim.obs$times[1:2])

obstx.mat = lapply(1:3, function(s) {
  dsdive.obstx.matrix(depth.bins = depth.bins, beta = beta, 
                      lambda = lambda, s0 = s, tstep = tstep, 
                      include.raw = TRUE, delta = 1e-10)
})

lambda.priors = list(
  c(4, 2 + 2/3), c(4, 13 + 1/3), c(4, 5)
)

beta.priors = list(
  c(31.5, 3.5), c(3.5, 31.5)
)

T1.prior.params = c(25, .04)
T2.prior.params = c(56, .06)

fit = dsdive.gibbs.obs(
  dsobs.list = list(sim.obs), t.stages.list = list(t.stages), 
  beta.init = beta, lambda.init = lambda, verbose = TRUE, maxit = 2, 
  pi1.prior = beta.priors[[1]], pi2.prior = beta.priors[[2]], 
  lambda1.prior = lambda.priors[[1]], lambda2.prior = lambda.priors[[2]], 
  lambda3.prior = lambda.priors[[3]], tstep = tstep, depth.bins = depth.bins, 
  T1.prior.params = T1.prior.params, T2.prior.params = T2.prior.params, 
  max.width = 100, max.width.offset = 30, t0.prior.params = c(1,1), 
  tf.prior.params = c(1,1), offsets = 0, offsets.tf = 0, warmup = 1, 
  delta = 1e-10)


detach(dive.sim$params)
detach(dive.sim)
}
