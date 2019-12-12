#' SMC estimation of likelihood for imputed dive trajectories
#'
#'
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
#' @param steps.max maximum number of transitions to sample before stopping, 
#'   regardless of whether \code{tf} is reached.
#' @param N number of SMC particles to maintain 
#' @param depths Dive bins in which the trajectory was observed
#' @param t Times at which depths were observed
#' @param tries.max maximum number of attempts to find a suitable ABC trajectory 
#'   before stopping sampling.  If \code{tries.max} is exceeded while sampling 
#'   any location, the returned log-likelihood estimate will be \code{-Inf}
#' @param eps ABC rejection sampling tolerance; maximum distance along graph 
#'   between simulated and observed states.  Set \code{eps=1} for exact 
#'   imputation.
#' @param dump.state if sampling fails and \code{dump.state==TRUE}, then the 
#'   particles in the particle filter will be returned, along with the timepoint 
#'   at which sampling failed
#' @param verbose if \code{TRUE}, then various diagnostic or informational 
#'   output will be printed while the ABC-SMC algorithm is running
#' @param n.samples the number of imputed trajectories to return.  If 
#'   \code{n.samples==0}, then no trajectories will be saved during sampling, 
#'   which can potentially yield faster inference because each step of the 
#'   ABC-SMC sampler requires fewer memory operations.
#' @param stage.init The dive stage in which all particles should be initialized
#' @param t0.dive Time at which dive started
#' @param shift.cfg vector: Optional arguments to bias sampling toward a 
#'   specific node.  See documentation for \code{dsdive.tx.params} for more 
#'   detail. c(depth.bias, rate.p).
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#'   
#' @example examples/ldabc.R
#' 
#' @export
#' 
#'
dsdive.ldabc = function(beta, lambda, sub.tx, surf.tx, depth.bins, 
                        steps.max = 1e3, N, depths, t, tries.max, eps,
                        dump.state = FALSE, verbose = FALSE, n.samples = 0,
                        stage.init = 1, t0.dive, shift.cfg = NULL, model) {
  
  # extract dimensional information
  nt = length(t)
  NN = N-1
  
  # initialize log-density with the constant term
  ld = (nt-1) * log(NN)
  ld.M = NULL
  
  #
  # initialize particle collections
  #
  
  particle = list(
    depths = depths[1],
    stages = stage.init,
    durations = NULL,
    times = t[1],
    logW = 0,
    t.stage2 = NA
  )
  class(particle) = 'dsdive'
  
  particles = list(
    resampling = rep(list(particle), NN),
    alive = list()
  )
  
  # sample particles for timepoint
  for(j in 2:nt) {
    
    if(verbose) {
      message(paste('Sampling timepoint', j))
    }
    
    # total sampling effort for timepoint
    M = 0
    
    # set biased sampling parameters
    if(is.null(shift.cfg)) {
      shift.params = NULL
    } else {
      shift.params = c(depths[j], shift.cfg[1], t[j], shift.cfg[2])
    }
    
    # compute sampling weights
    W = exp(scale(x = sapply(particles$resampling, function(p) p$logW), 
                  center = TRUE, scale = FALSE))
    
    # sample each particle
    for(i in 1:N) {
      
      err = eps
      
      # rejection sample to find acceptable particle
      for(m in 1:tries.max) {
        # resample an alive particle
        particle = particles$resampling[[sample(x = 1:NN, size = 1, prob = W)]]
        len = length(particle$depths)
        
        # propose particle; compute error
        p = dsdive.fwdsample(depth.bins = depth.bins, 
                             d0 = particle$depths[len], beta = beta, 
                             lambda = lambda, sub.tx = sub.tx, 
                             surf.tx = surf.tx, t0 = particle$times[len], 
                             tf = t[j], steps.max = steps.max, 
                             dur0 = particle$durations[len], 
                             s0 = particle$stages[len], t0.dive = t0.dive, 
                             shift.params = shift.params, 
                             t.stage2 = particle$t.stage2, model = model)
        
        # continue resampling if particle is not observable at required time
        if(length(p$times) > 0) {
          err = abs(p$depths[length(p$depths)] - depths[j])
        } else {
          err = eps
        }
        
        # stop sampling if acceptable particle is found
        if(err < eps)
          break
      }
      
      if(err < eps) {
        # update sampling effort for this timepoint
        M = M + m
        
        # save first N-1 particles
        if(i < N) {
          if(n.samples>0) { # update particle - save complete trajectory
            if(length(particle$durations) != length(particle$depths)) {
              particle$durations = c(particle$durations, p$durations)
            } else {
              particle$durations = c(particle$durations, p$durations[-1])
            }
            particle$depths = c(particle$depths, p$depths[-1])
            particle$times = c(particle$times, p$times[-1])
            particle$stages = c(particle$stages, p$stages[-1])
            
            
          } else { # update particle - only save last depth, duration, and stage
            particle$durations = p$durations[length(p$durations)]
            particle$depths = p$depths[length(p$depths)]
            particle$times = p$times[length(p$times)]
            particle$stages = p$stages[length(p$stages)]
          }
          
          # update stage 2 entry time, as necessary
          if(is.na(particle$t.stage2)) {
            stage2.inds = which(particles$stages==2)
            if(length(stage2.inds)>0) {
              particle$t.stage2 = particle$times[min(stage2.inds)]
            }
          }
          
          # save particle
          particle$logW = p$logW
          particles$alive[[i]] = particle
          
        }
      } else {
        # failed to find suitable proposal; fail sampling
        break
      }
    }
    
    ld.M = c(ld.M, M)
    
    if(err < eps) { 
      # update log-density estimate
      # lognbhdsize = log(2*eps)
      ld = ld - log(M-1) #- lognbhdsize
      
      # swap alive particles with resampling particles
      particles$resampling = particles$alive
      particles$alive = list()
    } else {
      # stop approximation if sampling has failed
      ld = -Inf
      break
    }
  }
  
  if(n.samples > NN) {
    warning('Requested number of samples exceeds size of saved particles;
             returning all saved particles instead.')
    n.samples = NN
  }
  
  # compute sampling weights
  W = exp(scale(x = sapply(particles$resampling, function(p) p$logW ), 
                center = TRUE, scale = FALSE))
  
  # package results
  res = list(
    ld = ld,
    sim = particles$resampling[sample(x = 1:NN, size = n.samples, 
                                      replace = FALSE, prob = W)],
    ld.M = ld.M
  )
  
  if(ld == -Inf) {
    if(dump.state) {
      res$particles = particles
      res$j = j
    }
  }
  res
}