# Simulate a bundle of Brownian Motion paths

Simulates `n_traj` independent brownian motion paths following the SDE:
\$\$dS_t = \mu dt + \sigma dW_t\$\$ With solution \$\$S_t = \mu t +
\sigma W_t\$\$ where \\W_t\\ is the Wiener process.

## Usage

``` r
brownian_matrix(dt, n_steps, n_traj, x0 = 0, mu = 0, sigma = 1)
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

[Wiener process](https://en.wikipedia.org/wiki/Wiener_process)

## See also

[`brownian_vector()`](https://antonikobosz.github.io/stochgen/reference/brownian_vector.md),
[`brownian_end_info()`](https://antonikobosz.github.io/stochgen/reference/brownian_end_info.md),

## Examples

``` r
# Simulate a bundle of Brownian Motion paths
M <- brownian_matrix(dt = 0.01, n_steps = 100, n_traj = 10)
matplot(M, type = "l", lty = 1,
        xlab = "step", ylab="S_t",
        main = "Brownian Motion")
```
