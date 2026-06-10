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
#' @param lambda intensity
#' @param n number of jumps
#' 
#' @description
#' Simulates a Poisson process with intensity \eqn{\lambda} by generating
#' the jump times \eqn{t_1, t_2, \ldots, t_n} of the first \code{n} jumps.
#' Inter-arrival times follow the Exponential distribution with parameter
#' \eqn{\lambda}.
#'
#' @references 
#' [Poisson point process](https://en.wikipedia.org/wiki/Poisson_point_process)
#' 
#' @examples 
#' # Simulate a Poisson process
#' lambda <- 2
#' n <- 50
#' times <- poiss_jump_times(lambda, n)
#' 
#' t_max <- max(times) + 1/lambda # extra padding
#' steps <- c(0, times, t_max)
#' N <- 0:length(times)  # process value
#' 
#' plot(steps, c(N, tail(N, 1)),
#'      type = "s",  # step function
#'      xlab = "t", ylab = "N(t)",
#'      main = expression("Poisson process, " ~ lambda == 2))
#' abline(a = 0, b = lambda, col = "red", lty = 2)
#' legend("topleft", legend = expression(lambda * t), col = "red", lty = 2)
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
#' where \eqn{W_t} is the Wiener process.
#'
#' @param dt      Time step size
#' @param n       Number of steps
#' @param x0      Initial value
#' @param mu      Drift coefficient
#' @param sigma   Volatility coefficient
#' @references
#' [Wiener process](https://en.wikipedia.org/wiki/Wiener_process)
#' @seealso
#' [brownian_matrix()], [brownian_end_info()],
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
#' where \eqn{W_t} is the Wiener process.
#'
#' @inheritParams brownian_vector
#'
#' @return
#' Numeric vector of length \code{n}
#' @examples
#' # Simulate a Geometric Brownian Motion path
#' path <- geom_brownian_vector(dt = 0.01, n = 100, x0 = 1, mu = 0, sigma = 0.2)
#' plot(path, type = "l")
#' @references 
#' [Geometric Brownian Motion](https://en.wikipedia.org/wiki/Geometric_Brownian_motion)
#' @seealso
#' [geom_brownian_matrix()], [geom_brownian_end_info()],
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
#' where \eqn{W_t} is the Wiener process.
#'
#' @inheritParams brownian_vector
#' @param theta   Reversion rate
#'
#' @return
#' Numeric vector of length \code{n}
#' @references 
#' [Ornstein–Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process)
#' @seealso
#' [ornstein_uhlenbeck_matrix()], [ornstein_uhlenbeck_end_info()],
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
#' where \eqn{W_t} is the Wiener process.
#'
#' @param dt      Time step size
#' @param n_steps Number of steps in each path
#' @param n_traj  Number of simulated paths
#' @param x0      Initial value
#' @param mu      Drift coefficient
#' @param sigma   Volatility coefficient
#' @references
#' [Wiener process](https://en.wikipedia.org/wiki/Wiener_process)
#' @seealso
#' [brownian_vector()], [brownian_end_info()],
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
#' where \eqn{W_t} is the Wiener process.
#'
#' @inheritParams brownian_matrix
#' @references 
#' [Geometric Brownian Motion](https://en.wikipedia.org/wiki/Geometric_Brownian_motion)
#' @seealso
#' [geom_brownian_vector()], [geom_brownian_end_info()],
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
#' where \eqn{W_t} is the Wiener process.
#'
#' @inheritParams brownian_matrix
#' @param theta   Reversion rate

#' @return
#' Numeric matrix of dimension \code{n_steps x n_traj} containing simulated paths.
#' Each column corresponds to one trajectory.
#' @references 
#' [Ornstein–Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process)
#' @seealso
#' [ornstein_uhlenbeck_vector()], [ornstein_uhlenbeck_end_info()],
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
#' Calculate final values of a bundle of brownian paths
#'
#' @description
#' Simulates \code{n_traj} independent brownian motion paths following the SDE:
#' \deqn{dS_t = \mu dt + \sigma dW_t}
#' With solution
#' \deqn{S_t = \mu t + \sigma W_t}
#' where \eqn{W_t} is the Wiener process.
#' For each path, returns its terminal value and whether it
#' exceeded the provided thresholds at any point.
#' 
#' @param dt      Time step size
#' @param n_steps Number of steps in each path
#' @param n_traj  Number of simulated paths
#' @param x0      Initial value
#' @param mu      Drift coefficient
#' @param sigma   Volatility coefficient
#' @param lower_thresh Lower threshold. A path is flagged if it falls below this value at any point
#' @param upper_thresh Upper threshold. A path is flagged if it exceeds this value at any point
#' 
#' @return A list with three elements:
#'   \describe{
#'     \item{xT}{A numeric vector of length \code{n_paths} containing the
#'       terminal value of each simulated path.}
#'     \item{crossed_upper}{A logical vector of length \code{n_paths}.
#'       \code{TRUE} if the path exceeded the upper threshold at any point.}
#'     \item{crossed_lower}{A logical vector of length \code{n_paths}.
#'       \code{TRUE} if the path fell below the lower threshold at any point.}
#'   }
#' 
#' @references
#' [Wiener process](https://en.wikipedia.org/wiki/Wiener_process)
#' @seealso
#' [brownian_vector()], [brownian_matrix()],
#' 
#' @examples
#' # Simulate 100 paths of standard Brownian motion (mu=0, sigma=1)
#' # over 100 steps of size 0.01 (total time T = 1)
#' result <- brownian_end_info(dt = 0.01, n_steps = 100, n_traj = 100)
#'
#' # Terminal values
#' hist(result$xT, main = "Distribution of S_T", xlab = "S_T")
#'
#' # Fraction of paths that crossed each threshold
#' mean(result$crossed_upper)
#' mean(result$crossed_lower)
#' @export
brownian_end_info <- function(dt, n_steps, n_traj, x0 = 0, mu = 0, sigma = 1,
                              lower_thresh = -1, upper_thresh = 1) {
  .Call("stochastic_end_info_r", x0, mu, sigma, 1, dt, lower_thresh,
        upper_thresh, n_steps, n_traj, .process_type$brownian,
        PACKAGE = "stochgen")
}

#' @title
#' Calculate final values of a bundle of geometric brownian paths
#' @description
#' Simulates \code{n_traj} independent geometric brownian motion paths following the SDE:
#' \deqn{dS_t = \mu S_t dt + \sigma S_t dW_t}
#' where \eqn{W_t} is the Wiener process.
#' For each path, returns its terminal value and whether it
#' exceeded the provided thresholds at any point.
#' 
#' @inheritParams brownian_end_info
#' 
#' @references 
#' [Geometric Brownian Motion](https://en.wikipedia.org/wiki/Geometric_Brownian_motion)
#' @seealso
#' [geom_brownian_vector()], [geom_brownian_matrix()],
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{xT}{A numeric vector of length \code{n_paths} containing the
#'       terminal value of each simulated path.}
#'     \item{crossed_upper}{A logical vector of length \code{n_paths}.
#'       \code{TRUE} if the path exceeded the upper threshold at any point.}
#'     \item{crossed_lower}{A logical vector of length \code{n_paths}.
#'       \code{TRUE} if the path fell below the lower threshold at any point.}
#'   }
#' 
#' @examples
#' # Simulate 100 paths of geometric Brownian motion (mu=0, sigma=1)
#' # over 100 steps of size 0.01 (total time T = 1)
#' result <- geom_brownian_end_info(
#'   dt = 0.01, n_steps = 100, n_traj = 100
#' )
#'
#' # Terminal values
#' hist(result$xT, main = "Distribution of S_T", xlab = "S_T")
#'
#' # Fraction of paths that crossed each threshold
#' mean(result$crossed_upper)
#' mean(result$crossed_lower)
#' @export
geom_brownian_end_info <- function(dt, n_steps, n_traj, x0 = 1, mu = 0,
                                   sigma = 1, lower_thresh = 0.5,
                                   upper_thresh = 2) {
  .Call("stochastic_end_info_r", x0, mu, sigma, 1, dt, lower_thresh,
        upper_thresh, n_steps, n_traj, .process_type$geom_brownian,
        PACKAGE = "stochgen")
}

#' @title
#' Calculate final values of a bundle ornstein-uhlenbeck processes 
#' @description
#' Simulates \code{n_traj} independent ornstein-uhlenbeck processes following the SDE:
#' \deqn{dS_t = -\theta(S_t-\mu)dt + \sigma dW_t}
#' where \eqn{W_t} is the Wiener process.
#' For each path, returns its terminal value and whether it
#' exceeded the provided thresholds at any point.
#' 
#' @inheritParams brownian_end_info
#' @param theta   Reversion rate
#' 
#' @references 
#' [Ornstein–Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process)
#' @seealso
#' [ornstein_uhlenbeck_vector()], [ornstein_uhlenbeck_matrix()],
#' @return A list with three elements:
#'   \describe{
#'     \item{xT}{A numeric vector of length \code{n_paths} containing the
#'       terminal value of each simulated path.}
#'     \item{crossed_upper}{A logical vector of length \code{n_paths}.
#'       \code{TRUE} if the path exceeded the upper threshold at any point.}
#'     \item{crossed_lower}{A logical vector of length \code{n_paths}.
#'       \code{TRUE} if the path fell below the lower threshold at any point.}
#'   }
#' @examples
#' # Simulate 100 paths of Ornstein-Uhlenbeck process (mu=0, theta=1, sigma=1)
#' # over 100 steps of size 0.01 (total time T = 1)
#' result <- ornstein_uhlenbeck_end_info(dt = 0.01, n_steps = 100, n_traj = 100)
#'
#' # Terminal values
#' hist(result$xT, main = "Distribution of S_T", xlab = "S_T")
#'
#' # Fraction of paths that crossed each threshold
#' mean(result$crossed_upper)
#' mean(result$crossed_lower)
#' 
#' @export
ornstein_uhlenbeck_end_info <- function(dt, n_steps, n_traj, x0 = 0,
                                        mu = 0, theta = 1, sigma = 1,
                                        lower_thresh = -1, upper_thresh = 1) {
  .Call("stochastic_end_info_r", x0, mu, sigma, theta, dt, lower_thresh,
        upper_thresh, n_steps, n_traj, .process_type$ornstein_uhlenbeck,
        PACKAGE = "stochgen")
}