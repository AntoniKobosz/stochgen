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
#' Simulate a Brownian Motion path
#'
#' @description
#' Simulates a brownian motion path following the SDE:
#' \deqn{dS_t = \mu dt + \sigma dW_t}
#' With solution
#' \deqn{S_t = \mu t + \sigma W_t}
#' where \eqn{W_t} is the Weiner process.
#'
#' @param dt      Time step size
#' @param n       Number of steps
#' @param x0      Initial value
#' @param mu      Drift coefficient
#' @param sigma   Volatility coefficient
#' @seealso
#' [brownian_matrix()], [brownian_end_info()],
#' [Wiener process](https://en.wikipedia.org/wiki/Wiener_process)
#'
#' @return
#' Numeric vector of length \code{n}
#' @examples
#' # Simulate a Brownian Motion path
#' path <- brownian_vector(dt = 0.01, n = 100, x0 = 0, mu = 0, sigma = 1)
#' plot(path, type = "l")
#' @export
brownian_vector <- function(dt, n, x0 = 0, mu = 0, sigma = 1) {
  .Call("stochastic_vector_r", x0, mu, sigma, 1, dt, n,
        .process_type$brownian, PACKAGE = "stochgen")
}
#' @title
#' Simulate a Geometric Brownian Motion path
#'
#' @description
#' Simulates a brownian motion path following the SDE:
#' \deqn{dS_t = \mu S_t dt + \sigma S_t dW_t}
#' where \eqn{W_t} is the Weiner process.
#'
#' @inheritParams brownian_vector
#'
#' @return
#' Numeric vector of length \code{n}
#' @examples
#' # Simulate a Geometric Brownian Motion path
#' path <- geom_brownian_vector(dt = 0.01, n = 100, x0 = 1, mu = 0, sigma = 0.2)
#' plot(path, type = "l")
#' @seealso
#' [geom_brownian_matrix()], [geom_brownian_end_info()],
#' [Geometric Brownian Motion](https://en.wikipedia.org/wiki/Geometric_Brownian_motion)
#'
#' @export
geom_brownian_vector <- function(dt, n, x0 = 1, mu = 0, sigma = 1) {
  .Call("stochastic_vector_r", x0, mu, sigma, 1, dt, n,
        .process_type$geom_brownian, PACKAGE = "stochgen")
}

#' @title
#' Simulate an Ornstein-Uhlenbeck path
#'
#' @description
#' Simulates an Ornstein-Uhlenbeck path following the SDE:
#' \deqn{dS_t = -\theta(S_t-\mu)dt + \sigma dW_t}
#' where \eqn{W_t} is the Weiner process.
#'
#' @inheritParams brownian_vector
#' @param theta   Reversion rate
#'
#' @return
#' Numeric vector of length \code{n}
#' @seealso
#' [ornstein_uhlenbeck_matrix()], [ornstein_uhlenbeck_end_info()],
#' [Ornstein–Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process)
#'
#' @examples
#' # Simulate an Ornstein-Uhlenbeck path
#' path <- ornstein_uhlenbeck_vector(dt = 0.01, n = 100, x0 = 0,
#'                                   mu = 0, theta = 1, sigma = 1)
#' plot(path, type = "l")

#' @export
ornstein_uhlenbeck_vector <- function(dt, n, x0 = 0, mu = 0,
                                      theta = 1, sigma = 1) {
  .Call("stochastic_vector_r", x0, mu, sigma, theta, dt, n,
        .process_type$ornstein_uhlenbeck, PACKAGE = "stochgen")
}




#' @title
#' Simulate a bundle of Brownian Motion paths
#'
#' @description
#' Simulates \code{n_traj} independent brownian motion paths following the SDE:
#' \deqn{dS_t = \mu dt + \sigma dW_t}
#' With solution
#' \deqn{S_t = \mu t + \sigma W_t}
#' where \eqn{W_t} is the Weiner process.
#'
#' @param dt      Time step size
#' @param n_steps Number of steps in each path
#' @param n_traj  Number of simulated paths
#' @param x0      Initial value
#' @param mu      Drift coefficient
#' @param sigma   Volatility coefficient
#' @seealso
#' [brownian_vector()], [brownian_end_info()],
#' [Wiener process](https://en.wikipedia.org/wiki/Wiener_process)
#'
#' @return
#' Numeric matrix of dimension \code{n_steps x n_traj} containing simulated paths.
#' Each column corresponds to one trajectory.
#' @examples
#' # Simulate a bundle of Brownian Motion paths
#' M <- brownian_matrix(dt = 0.01, n_steps = 100, n_traj = 10)
#' matplot(M, type = "l", lty = 1,
#'         xlab = "step", ylab="S_t",
#'         main = "Brownian Motion")
#' @export
brownian_matrix <- function(dt, n_steps, n_traj, x0 = 0, mu = 0, sigma = 1) {
  .Call("stochastic_matrix_r", x0, mu, sigma, 1, dt, n_steps, n_traj,
        .process_type$brownian, PACKAGE = "stochgen")
}

#' @title
#' Simulate a bundle of Geometric Brownian Motion paths
#'
#' @description
#' Simulates \code{n_traj} independent geometric brownian motion paths following the SDE:
#' \deqn{dS_t = \mu S_t dt + \sigma S_t dW_t}
#' where \eqn{W_t} is the Weiner process.
#'
#' @inheritParams brownian_matrix
#' @seealso
#' [geom_brownian_vector()], [geom_brownian_end_info()],
#' [Geometric Brownian Motion](https://en.wikipedia.org/wiki/Geometric_Brownian_motion)
#'
#' @return
#' Numeric matrix of dimension \code{n_steps x n_traj} containing simulated paths.
#' Each column corresponds to one trajectory.
#' @examples
#' # Simulate a bundle of Geometric Brownian Motion paths
#' M <- geom_brownian_matrix(dt = 0.01, n_steps = 100, n_traj = 10)
#' matplot(M, type = "l", lty = 1,
#'         xlab = "step", ylab="S_t",
#'         main = "Geometric Brownian Motion")
#' @export
geom_brownian_matrix <- function(dt, n_steps, n_traj, x0 = 1,
                                 mu = 0, sigma = 1) {
  .Call("stochastic_matrix_r", x0, mu, sigma, 1, dt, n_steps, n_traj,
        .process_type$geom_brownian, PACKAGE = "stochgen")
}

#' @title
#' Simulate a bundle of Ornstein-Uhlenbeck paths
#'
#' @description
#' Simulates \code{n_traj} independent Ornstein-Uhlenbeck paths following the SDE:
#' \deqn{dS_t = -\theta(S_t-\mu)dt + \sigma dW_t}
#' where \eqn{W_t} is the Weiner process.
#'
#' @inheritParams brownian_matrix
#' @param theta   Reversion rate

#' @return
#' Numeric matrix of dimension \code{n_steps x n_traj} containing simulated paths.
#' Each column corresponds to one trajectory.
#' @seealso
#' [ornstein_uhlenbeck_vector()], [ornstein_uhlenbeck_end_info()],
#' [Ornstein–Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process)
#'
#' @examples
#' # Simulate a bundle of Ornstein-Uhlenbeck paths
#' M <- ornstein_uhlenbeck_matrix(dt = 0.01, n_steps = 100, n_traj = 10)
#' matplot(M, type = "l", lty = 1,
#'         xlab = "step", ylab="S_t",
#'         main = "Ornstein-Uhlenbeck paths")
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