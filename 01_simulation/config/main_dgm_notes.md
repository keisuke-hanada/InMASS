# Main Simulation DGM Notes

These notes document the stage-one v2 implementation for the current main scenarios.

## External Sample Sizes

For each scenario, v2 generates one external total sample size per external study as:

```r
n_external <- 2L * round(runif(K, n / 2, 2 * n))
```

This matches the legacy C++ implementation. The resulting external study totals are even integers in `[n, 4n]`, so each external randomized study has balanced treatment and control arms with `n_external / 2` participants per arm. This is consistent with the manuscript statement that external trial sample sizes are generated from `U(n, 4n)`, with the implementation detail that totals are rounded to even integers for 1:1 external arm allocation.

## Covariate Distributions

The normal covariate scenario uses:

```r
X = Z + mu_k,  Z ~ N(0, 1)
```

so `X ~ N(mu_k, 1)`.

The chi-squared covariate scenario uses the legacy standardized two-normal construction:

```r
X = (Z1^2 + Z2^2) / 2 + mu_k - 1,  Z1, Z2 ~ N(0, 1)
```

Since `(Z1^2 + Z2^2) / 2` has mean 1 and variance 1, this gives mean `mu_k` and variance 1. This matches the current main-simulation chi-squared covariate DGM used by v1.

## Outcome Model

The current main one-covariate outcome model is:

```r
Y = 1 + 2 * Z - X + 0.5 * Z * X + error
```

with `error ~ N(0, sigma^2)` and target covariate mean zero. Therefore the target-population average treatment effect is `true_delta = 2`.
