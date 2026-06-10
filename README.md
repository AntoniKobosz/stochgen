# stochgen

Simulation of stochastic processes in R.

The package provides lightweight tools for simulating common stochastic processes such as:
Brownian motion, Geometric Brownian Motion, Ornstein–Uhlenbeck process, and Poisson jump processes.

Core computations are implemented in C with GSL and OpenMP for performance.

---

For more information, visit:
https://antonikobosz.github.io/stochgen/

---
## Installation

### From GitHub (development version)

```r
install.packages("remotes")

remotes::install_github("AntoniKobosz/stochgen")
```

---
# License
MIT, see LICENSE.