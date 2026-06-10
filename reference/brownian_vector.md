# Simulate a Brownian Motion path

Simulates a brownian motion path following the SDE: \$\$dS_t = \mu dt +
\sigma dW_t\$\$ With solution \$\$S_t = \mu t + \sigma W_t\$\$ where
\\W_t\\ is the Wiener process.

## Usage

``` r
brownian_vector(dt, n, x0 = 0, mu = 0, sigma = 1)
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

[Wiener process](https://en.wikipedia.org/wiki/Wiener_process)

## See also

[`brownian_matrix()`](https://antonikobosz.github.io/stochgen/reference/brownian_matrix.md),
[`brownian_end_info()`](https://antonikobosz.github.io/stochgen/reference/brownian_end_info.md),

## Examples

``` r
# Simulate a Brownian Motion path
path <- brownian_vector(dt = 0.01, n = 100, x0 = 0, mu = 0, sigma = 1)
plot(path, type = "l")
```
