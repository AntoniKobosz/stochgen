# Simulate an Ornstein-Uhlenbeck path

Simulates an Ornstein-Uhlenbeck path following the SDE: \$\$dS_t =
-\theta(S_t-\mu)dt + \sigma dW_t\$\$ where \\W_t\\ is the Wiener
process.

## Usage

``` r
ornstein_uhlenbeck_vector(dt, n, x0 = 0, mu = 0, theta = 1, sigma = 1)
```

## Arguments

- dt:

  Time step size

- n:

  Number of steps

- x0:

  Initial value

- mu:

  Drift coefficient

- theta:

  Reversion rate

- sigma:

  Volatility coefficient

## Value

Numeric vector of length `n`

## References

[Ornstein–Uhlenbeck
process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process)

## See also

[`ornstein_uhlenbeck_matrix()`](https://antonikobosz.github.io/stochgen/reference/ornstein_uhlenbeck_matrix.md),
[`ornstein_uhlenbeck_end_info()`](https://antonikobosz.github.io/stochgen/reference/ornstein_uhlenbeck_end_info.md),

## Examples

``` r
# Simulate an Ornstein-Uhlenbeck path
path <- ornstein_uhlenbeck_vector(dt = 0.01, n = 100, x0 = 0,
                                  mu = 0, theta = 1, sigma = 1)
plot(path, type = "l")
```
