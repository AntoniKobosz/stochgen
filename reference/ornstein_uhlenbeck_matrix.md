# Simulate a bundle of Ornstein-Uhlenbeck paths

Simulates `n_traj` independent Ornstein-Uhlenbeck paths following the
SDE: \$\$dS_t = -\theta(S_t-\mu)dt + \sigma dW_t\$\$ where \\W_t\\ is
the Wiener process.

## Usage

``` r
ornstein_uhlenbeck_matrix(
  dt,
  n_steps,
  n_traj,
  x0 = 0,
  mu = 0,
  theta = 1,
  sigma = 1
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

## Value

Numeric matrix of dimension `n_steps x n_traj` containing simulated
paths. Each column corresponds to one trajectory.

## References

[Ornstein–Uhlenbeck
process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process)

## See also

[`ornstein_uhlenbeck_vector()`](https://antonikobosz.github.io/stochgen/reference/ornstein_uhlenbeck_vector.md),
[`ornstein_uhlenbeck_end_info()`](https://antonikobosz.github.io/stochgen/reference/ornstein_uhlenbeck_end_info.md),

## Examples

``` r
# Simulate a bundle of Ornstein-Uhlenbeck paths
M <- ornstein_uhlenbeck_matrix(dt = 0.01, n_steps = 100, n_traj = 10)
matplot(M, type = "l", lty = 1,
        xlab = "step", ylab="S_t",
        main = "Ornstein-Uhlenbeck paths")
```
