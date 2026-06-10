# Calculate final values of a bundle of geometric brownian paths

Simulates `n_traj` independent geometric brownian motion paths following
the SDE: \$\$dS_t = \mu S_t dt + \sigma S_t dW_t\$\$ where \\W_t\\ is
the Wiener process. For each path, returns its terminal value and
whether it exceeded the provided thresholds at any point.

## Usage

``` r
geom_brownian_end_info(
  dt,
  n_steps,
  n_traj,
  x0 = 1,
  mu = 0,
  sigma = 1,
  lower_thresh = 0.5,
  upper_thresh = 2
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

[Geometric Brownian
Motion](https://en.wikipedia.org/wiki/Geometric_Brownian_motion)

## See also

[`geom_brownian_vector()`](https://antonikobosz.github.io/stochgen/reference/geom_brownian_vector.md),
[`geom_brownian_matrix()`](https://antonikobosz.github.io/stochgen/reference/geom_brownian_matrix.md),

## Examples

``` r
# Simulate 100 paths of geometric Brownian motion (mu=0, sigma=1)
# over 100 steps of size 0.01 (total time T = 1)
result <- geom_brownian_end_info(
  dt = 0.01, n_steps = 100, n_traj = 100
)

# Terminal values
hist(result$xT, main = "Distribution of S_T", xlab = "S_T")


# Fraction of paths that crossed each threshold
mean(result$crossed_upper)
#> [1] 0.26
mean(result$crossed_lower)
#> [1] 0.67
```
