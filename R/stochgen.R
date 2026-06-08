#' @useDynLib stochgen, .registration = TRUE
"_PACKAGE"


#' @title
#' Generete a vector of jump times of a Poisson process
#'
#' @description
#' A Poisson process with intensity \eqn{\lambda} can be constructed by
#' counting jumps which occur at random such that the time between jumps
#' follows the Exponential distribution with parameter \eqn{\lambda}
#'
#' @export
poiss_jump_times <- function(lambda, n) {
  .Call("poiss_jump_times_r", lambda, n, PACKAGE = "stochgen")
}

#' @title
#' TODO
#' 
#' @description
#' TODO
#' 
#' @param
#' TODO
#' 
#' @seealso
#' TODO
#' 
#' @export
brownian_vector <- function(x0, mu, sigma, dt, n) {
  .Call("brownian_vector_r", x0, mu, sigma, dt, n, PACKAGE = "stochgen")
}

#' @title
#' TODO
#' 
#' @description
#' TODO
#' 
#' @param
#' TODO
#' 
#' @seealso
#' TODO
#' 
#' @export
brownian_matrix <- function(x0, mu, sigma, dt, n, m) {
  .Call("brownian_matrix_r", x0, mu, sigma, dt, n, m, PACKAGE = "stochgen")
}

#' @title
#' TODO
#' 
#' @description
#' TODO
#' 
#' @param
#' TODO
#' 
#' @seealso
#' TODO
#' 
#' @export
brownian_end_info <- function(x0, mu, sigma, dt,
                              lower_thresh, upper_thresh, n, m) {
  .Call("brownian_end_info_r", x0, mu, sigma, dt,
        lower_thresh, upper_thresh, n, m, PACKAGE = "stochgen")
}