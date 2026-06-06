#' @useDynLib stochgen, .registration = TRUE
"_PACKAGE"


make_seed <- function() {
    sample(.Machine$integer.max, 1)
}


#' @title 
#' Generete a vector of jump times of a Poisson process
#' 
#' @description
#' A Poisson process with intensity \eqn{\lambda} can be constructed by counting jumps
#' which occur at random such that the time between jumps
#' follows the Exponential distribution with parameter \lambda
#' 
#' @export 
poiss_jump_times <- function(lambda, n) {
    seed <- make_seed()
    .Call("poiss_jump_times_r",lambda,n,seed,PACKAGE="stochgen")
}
