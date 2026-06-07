#' @useDynLib stochgen, .registration = TRUE
"_PACKAGE"


#' @title
#' Generete a vector of jump times of a Poisson process
#'
#' @description
#' A Poisson process with intensity \eqn{\lambda} can be constructed by
#' counting jumps which occur at random such that the time between jumps
#' follows the Exponential distribution with parameter \lambda
#'
#' @export
poiss_jump_times <- function(lambda, n) {
  .Call("poiss_jump_times_r", lambda, n, PACKAGE = "stochgen")
}
