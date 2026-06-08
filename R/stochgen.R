#' @useDynLib stochgen, .registration = TRUE
"_PACKAGE"


.process_type <- list(
  brownian = 0L,
  geom_brownian = 1L,
  ornstein_uhlenbeck = 2L
)

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
brownian_vector <- function(dt, n, x0 = 0, mu = 0, sigma = 1) {
  .Call("stochastic_vector_r", x0, mu, sigma, 1, dt, n,
        .process_type$brownian, PACKAGE = "stochgen")
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
geom_brownian_vector <- function(dt, n, x0 = 1, mu = 0, sigma = 1) {
  .Call("stochastic_vector_r", x0, mu, sigma, 1, dt, n,
        .process_type$geom_brownian, PACKAGE = "stochgen")
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
ornstein_uhlenbeck_vector <- function(dt, n, x0 = 0, mu = 0,
                                      theta = 1, sigma = 1) {
  .Call("stochastic_vector_r", x0, mu, sigma, theta, dt, n,
        .process_type$ornstein_uhlenbeck, PACKAGE = "stochgen")
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
brownian_matrix <- function(dt, n_steps, n_traj, x0 = 0, mu = 0, sigma = 1) {
  .Call("stochastic_matrix_r", x0, mu, sigma, 1, dt, n_steps, n_traj,
        .process_type$brownian, PACKAGE = "stochgen")
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
geom_brownian_matrix <- function(dt, n_steps, n_traj, x0 = 1,
                                 mu = 0, sigma = 1) {
  .Call("stochastic_matrix_r", x0, mu, sigma, 1, dt, n_steps, n_traj,
        .process_type$geom_brownian, PACKAGE = "stochgen")
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
ornstein_uhlenbeck_matrix <- function(dt, n_steps, n_traj, x0 = 0,
                                      mu = 0, theta = 1, sigma = 1) {
  .Call("stochastic_matrix_r", x0, mu, sigma, theta, dt, n_steps, n_traj,
        .process_type$ornstein_uhlenbeck, PACKAGE = "stochgen")
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
brownian_end_info <- function(dt, n_steps, n_traj, x0 = 0, mu = 0, sigma = 1,
                              lower_thresh = -1, upper_thresh = 1) {
  .Call("stochastic_end_info_r", x0, mu, sigma, 1, dt, lower_thresh,
        upper_thresh, n_steps, n_traj, .process_type$brownian,
        PACKAGE = "stochgen")
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
geom_brownian_end_info <- function(dt, n_steps, n_traj, x0 = 1, mu = 0,
                                   sigma = 1, lower_thresh, upper_thresh) {
  .Call("stochastic_end_info_r", x0, mu, sigma, 1, dt, lower_thresh,
        upper_thresh, n_steps, n_traj, .process_type$geom_brownian,
        PACKAGE = "stochgen")
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
ornstein_uhlenbeck_end_info <- function(dt, n_steps, n_traj, x0 = 0,
                                        mu = 0, theta = 1, sigma = 1,
                                        lower_thresh = -1, upper_thresh = 1) {
  .Call("stochastic_end_info_r", x0, mu, sigma, theta, dt, lower_thresh,
        upper_thresh, n_steps, n_traj, .process_type$ornstein_uhlenbeck,
        PACKAGE = "stochgen")
}