#' The geometry of humans' knowledge of navigation space.
#'
#' A dataset from a research by Warren et.al. (2017) on the geometry of humans'
#' knowledge of navigation space.
#'
#' In their research Warren et.al. (2017) conduct an experiment in which a total
#' of 20 participants had to navigate through one of two versions of a virtual
#' maze using virtual reality headsets. One is the standard version, Euclidean
#' maze, with several objects hidden inside and shortcuts to these objects. The
#' other maze is exactly the same apart from that it contains wormholes by which
#' participants can be 'teleported' from one place in the maze to another, the
#' non-Euclidean maze.
#'
#' The task participants had to execute was to walk to the several target
#' objects in the maze from a certain start object. There were two types of
#' targets, probe and standard and for each of the types there were 4 different
#' ones resulting in a total of 8 targets.
#'
#' Before the actual experiment participants first completed a training phase in
#' which participants were trained to walk to each target object.
#'
#' The variable of interest in this experiment was the angular error, the angle
#' between the intial walking direction of the participant and the target's
#' location (for specifics see article by Warren et.al.).
#'
#' The angular error is a circular variable and can be described and analyzed
#' using circular statistics.
#'
#' @format A data frame with 160 rows and 7 variables:
#' \describe{
#'   \item{Subject}{a numeric variable indicating the participant number}
#'   \item{Trial.no}{a numeric variable indicating the trial number of the
#'   target a participant had to locate (1-8)}
#'   \item{Maze}{a between subjects factor variable indicating the type of maze
#'   a participant was in (0 = Euclidean, 1 = non-Euclidean)}
#'   \item{Trial.type}{a within subjects factor indicating the type of target (0
#'   = standard, 1 = probe)}
#'    \item{Error}{a numeric variable containing the angular error in degrees}
#'   \item{Learn}{a numeric variable indicating the number of trials a
#'   participant completed in the training phase}
#'   \item{Learn.c}{mean centered Learning}
#'   \item{Error.rad}{a numeric variable containing the angular error in radians}
#'}
#'
#' @source \url{https://doi.org/10.1016/j.cognition.2017.05.020}
#'
"Maps"

#' Phase differences in hand flexion-extension movements.
#'
#' A dataset from a research by Puglisi et.al. (2017) on the role of attention
#' in human motor resonance.
#'
#' In their research Puglisi et.al. (2017) conduct a between subjects experiment
#' in which `observers' in multiple degrees of explicitness are asked to look at
#' the movement of a hand of the `mover' or other object in order to evaluate
#' the role of attention in motor resonant response.
#'
#' The experiment has four conditions:
#'
#' 1. The `explicit observation' condition (n = 14), where observers are
#' explicitly instructed to observe the hand.
#'
#' 2. The `semi-implicit observation' condition (n = 14) where the observers
#' have to do a task that requires implicit observation of the hand.
#'
#' 3. The `implicit observation' condition (n = 14), where observers have to do
#' a task that is independent of the observation of the hand that is moving in
#' front of them.
#'
#' 4. A baseline condition (n=14) where there is no moving hand but observers
#' have to look at an inanimate object that moves in an identical manner to the
#' hand.
#'
#' The idea of motor resonance is then that the `observer' because he is looking
#' explicitly or implicitly at the hand of the `mover' starts moving his or her
#' hand in the same manner. This is the resonant response. In each condition the
#' hand movements of the observers were measured and the phase difference
#' between the observers' hand and the hand they observed (the movers' hand) was
#' calculated. This was not done for the baseline condition because in this
#' condition there was no periodic pattern of movement in the `observers' hand.
#'
#' The phase difference is a circular variable and can be described and analyzed
#' using circular statistics.
#'
#' @format A data frame with 42 rows and 7 variables:
#' \describe{
#'   \item{Cond}{a factor variable indicating the condition a participant
#'   was placed in; 1 = 'implicit', 2 = 'semi.implicit' or 3 = 'explicit'}
#'   \item{PhaseDiff}{a numeric variable the phase difference between 'observer'
#'   and 'mover' in degrees}
#'   \item{AvAmp}{a numeric variable indicating the average amplitude of the
#'   hand movement of the 'observer'}
#'   \item{Phaserad}{a numeric variable the phase difference between 'observer'
#'   and 'mover' in radians}
#' }
#'
#' @source \url{https://doi.org/10.1371/journal.pone.0177457}
#'
"Motor"


