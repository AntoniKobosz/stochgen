# Simulate a bundle of Geometric Brownian Motion paths

Simulates `n_traj` independent geometric brownian motion paths following
the SDE: \$\$dS_t = \mu S_t dt + \sigma S_t dW_t\$\$ where \\W_t\\ is
the Wiener process.

## Usage

``` r
geom_brownian_matrix(dt, n_steps, n_traj, x0 = 1, mu = 0, sigma = 1)
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

## Value

Numeric matrix of dimension `n_steps x n_traj` containing simulated
paths. Each column corresponds to one trajectory.

## References

[Geometric Brownian
Motion](https://en.wikipedia.org/wiki/Geometric_Brownian_motion)

## See also

[`geom_brownian_vector()`](https://antonikobosz.github.io/stochgen/reference/geom_brownian_vector.md),
[`geom_brownian_end_info()`](https://antonikobosz.github.io/stochgen/reference/geom_brownian_end_info.md),

## Examples

``` r
# Simulate a bundle of Geometric Brownian Motion paths
M <- geom_brownian_matrix(dt = 0.01, n_steps = 100, n_traj = 10)
matplot(M, type = "l", lty = 1,
        xlab = "step", ylab="S_t",
        main = "Geometric Brownian Motion")
```
