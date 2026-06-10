# Simulate a Geometric Brownian Motion path

Simulates a brownian motion path following the SDE: \$\$dS_t = \mu S_t
dt + \sigma S_t dW_t\$\$ where \\W_t\\ is the Wiener process.

## Usage

``` r
geom_brownian_vector(dt, n, x0 = 1, mu = 0, sigma = 1)
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

- sigma:

  Volatility coefficient

## Value

Numeric vector of length `n`

## References

[Geometric Brownian
Motion](https://en.wikipedia.org/wiki/Geometric_Brownian_motion)

## See also

[`geom_brownian_matrix()`](https://antonikobosz.github.io/stochgen/reference/geom_brownian_matrix.md),
[`geom_brownian_end_info()`](https://antonikobosz.github.io/stochgen/reference/geom_brownian_end_info.md),

## Examples

``` r
# Simulate a Geometric Brownian Motion path
path <- geom_brownian_vector(dt = 0.01, n = 100, x0 = 1, mu = 0, sigma = 0.2)
plot(path, type = "l")
```
