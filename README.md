# stochgen

Simulation of stochastic processes in R.

The package provides lightweight tools for simulating common stochastic processes such as Brownian motion,
Geometric brownian motion, Ornstein-Uhlenbeck, and Poisson jump processes. Core computations are implemented in C with  use of GSL and OpenMP for performance.
\\
For more information, check out the 
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