# Multi-Covariate Robustness DGM Notes

This simulation family is separate from the original main scenarios and is labeled:

```r
simulation_family = "robustness_multicov"
dgm = "multicov"
```

The outcome model is:

```r
Y = beta0 + delta_T * z + beta1 * x1 + beta2 * x2 +
    beta3 * z * x1 + beta4 * z * x2 + error
```

with `error ~ N(0, 1)` and:

```r
(delta_T, beta0, beta1, beta2, beta3, beta4) =
  (2, 1, -1, 0.5, 0.5, -0.25)
```

External study covariates are generated as:

```r
x1 ~ N(mu_1k, 1),  mu_1k = 4 * (k - 1) / (K - 1) - 1
x2 ~ N(mu_2k, 1),  mu_2k = 2 * (k - 1) / (K - 1) - 0.5
```

The target trial uses `mu_1T = 0` and `mu_2T = 0`, so the target-population average treatment effect is `true_delta = 2`.

The correctly specified working model is:

```r
yik ~ 1 + z + x1 + x2 + z:x1 + z:x2
```

The misspecified working model is:

```r
yik ~ 1 + z
```

The current global InMASS density-ratio implementation uses the reconstructed covariates and their squared terms in the propensity-score model. For this DGM, diagnostics record:

```r
density_ratio_covariates = "x1+x2"
density_ratio_quadratic_terms = "I(x1^2)+I(x2^2)"
```
