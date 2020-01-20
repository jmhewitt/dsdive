#' Compute transition matrix for Markov chain representation of dive model
#' 
#' Given model parameters and transition times, the 3 stage dive model can be 
#' represented as a Markov chain over a state space where each state records
#' the last depth bin, current depth bin, and current dive stage.  This 
#' function computes the complete transition matrix for the state space, given 
#' model parameters and time at which the transition will take place.
#'   
#' @param t0 time at which transition parameters should be computed
#' @param depth.bins \eqn{n x 2} Matrix that defines the depth bins.  The first 
#'   column defines the depth at the center of each depth bin, and the second 
#'   column defines the half-width of each bin.
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
#' @param inflation.factor.lambda In order to facilitate bridged transitions, 
#'   the transition rate of the overall process must be inflated to allow the 
#'   possibility of self-transitions.  Self-transitions allow bridged paths to 
#'   dynamically modify the total number of transitions between observed values
#'   so that a valid path between observations is always possible.  The 
#'   \code{inflation.factor.lambda} parameter implicitly controls the number of 
#'   self-transitions that will occur.  Larger values will create more 
#'   self-transitions.
#' @param min.depth As a computational efficiency option, only compute 
#'   transition parameters for depth bins at and above \code{min.depth}.
#' @param max.depth As a computational efficiency option, only compute 
#'   transition parameters for depth bins at and below \code{max.depth}.
#' @param t0.dive Time at which dive started
#' @param lambda.max Arrival rate for the parent Poisson process that will
#'   be thinned.  \code{lambda.max} will be scaled by 
#' @param t.stage2 time at which second stage was entered
#' @param model Either \code{"conditional"} or \code{"logit"} depending on the 
#'   method used to determine stage transition probability curves
#'   
#' @example examples/min.tx.R
#' 
#' @param d0 starting index
#' @param df ending index
#' @param ns number of stage transitions
#' 
#' @return minimum number of up/down transitions needed to go from d0 to df 
#'  allowing for ns stage transitions
#' 
#' @export
#' 
min.tx = function(d0, df, ns) {
  # distance to travel
 dtx = abs(df-d0) 
 # determine number of transitions needed
 if(dtx >= ns) {
   dtx
 } else {
   # determine even/odd, counting 0 as even
   dtx.even = dtx %% 2 == 0
   ns.even = ns %% 2 == 0
   # require extra transitions to satisfy parity requirements
   ns + ifelse(dts.even == ns.even, 0, 1)
 }
}