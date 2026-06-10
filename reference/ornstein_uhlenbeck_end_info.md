# Calculate final values of a bundle ornstein-uhlenbeck processes

Simulates `n_traj` independent ornstein-uhlenbeck processes following
the SDE: \$\$dS_t = -\theta(S_t-\mu)dt + \sigma dW_t\$\$ where \\W_t\\
is the Wiener process. For each path, returns its terminal value and
whether it exceeded the provided thresholds at any point.

## Usage

``` r
ornstein_uhlenbeck_end_info(
  dt,
  n_steps,
  n_traj,
  x0 = 0,
  mu = 0,
  theta = 1,
  sigma = 1,
  lower_thresh = -1,
  upper_thresh = 1
)
```

## Arguments

- dt:

  Time step size

- n_steps:

  Number of steps in each path

- n_traj:

  Number of simulated paths

- x0:

  Initial value

- mu:

  Drift coefficient

- theta:

  Reversion rate

- sigma:

  Volatility coefficient

- lower_thresh:

  Lower threshold. A path is flagged if it falls below this value at any
  point

- upper_thresh:

  Upper threshold. A path is flagged if it exceeds this value at any
  point

## Value

A list with three elements:

- xT:

  A numeric vector of length `n_paths` containing the terminal value of
  each simulated path.

- crossed_upper:

  A logical vector of length `n_paths`. `TRUE` if the path exceeded the
  upper threshold at any point.

- crossed_lower:

  A logical vector of length `n_paths`. `TRUE` if the path fell below
  the lower threshold at any point.

## References

[Ornstein–Uhlenbeck
process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process)

## See also

[`ornstein_uhlenbeck_vector()`](https://antonikobosz.github.io/stochgen/reference/ornstein_uhlenbeck_vector.md),
[`ornstein_uhlenbeck_matrix()`](https://antonikobosz.github.io/stochgen/reference/ornstein_uhlenbeck_matrix.md),

## Examples

``` r
# Simulate 100 paths of Ornstein-Uhlenbeck process (mu=0, theta=1, sigma=1)
# over 100 steps of size 0.01 (total time T = 1)
result <- ornstein_uhlenbeck_end_info(dt = 0.01, n_steps = 100, n_traj = 100)

# Terminal values
hist(result$xT, main = "Distribution of S_T", xlab = "S_T")


# Fraction of paths that crossed each threshold
mean(result$crossed_upper)
#> [1] 0.11
mean(result$crossed_lower)
#> [1] 0.21
```
