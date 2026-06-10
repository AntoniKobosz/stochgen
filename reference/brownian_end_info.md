# Calculate final values of a bundle of brownian paths

Simulates `n_traj` independent brownian motion paths following the SDE:
\$\$dS_t = \mu dt + \sigma dW_t\$\$ With solution \$\$S_t = \mu t +
\sigma W_t\$\$ where \\W_t\\ is the Wiener process. For each path,
returns its terminal value and whether it exceeded the provided
thresholds at any point.

## Usage

``` r
brownian_end_info(
  dt,
  n_steps,
  n_traj,
  x0 = 0,
  mu = 0,
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

[Wiener process](https://en.wikipedia.org/wiki/Wiener_process)

## See also

[`brownian_vector()`](https://antonikobosz.github.io/stochgen/reference/brownian_vector.md),
[`brownian_matrix()`](https://antonikobosz.github.io/stochgen/reference/brownian_matrix.md),

## Examples

``` r
# Simulate 100 paths of standard Brownian motion (mu=0, sigma=1)
# over 100 steps of size 0.01 (total time T = 1)
result <- brownian_end_info(dt = 0.01, n_steps = 100, n_traj = 100)

# Terminal values
hist(result$xT, main = "Distribution of S_T", xlab = "S_T")


# Fraction of paths that crossed each threshold
mean(result$crossed_upper)
#> [1] 0.3
mean(result$crossed_lower)
#> [1] 0.31
```
