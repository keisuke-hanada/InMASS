# Nonlinear Treatment-Effect Heterogeneity Scenario

This reviewer-requested robustness scenario is kept separate from the main simulation and from the multi-covariate robustness scenario.

simulation_family = "robustness_nonlinear"
dgm = "nonlinear"

## Data-Generating Model

For participant i in external trial S_k, k = 1, ..., K:

Y_ki = beta0 + delta_T z_ki + beta1 x_ki + beta2 x_ki^2
       + beta3 z_ki x_ki + beta4 z_ki (x_ki^2 - 1) + epsilon_ki

epsilon_ki ~ N(0, 1)

x_ki ~ N(mu_k, 1)

mu_k = 4(k - 1) / (K - 1) - 1

The target trial uses mu_T = 0.

Parameter values:

(delta_T, beta0, beta1, beta2, beta3, beta4) = (2, 1, -1, 0.25, 0.5, 0.25)

Because E_T[X] = 0 and E_T[X^2 - 1] = 0 under the target covariate distribution, true_delta = 2.

## Aggregate Data

Aggregate data include arm-specific sample size, outcome mean and variance, and x mean and variance.

The nonlinear aggregate features used in the analysis are derived from the reported mean and variance:

E[X^2] = Var(X) + E[X]^2

E[X^2 - 1] = Var(X) + E[X]^2 - 1

The analysis path does not use hidden external IPD to compute these nonlinear aggregate moments.

## Working Models

Correctly specified:

yik ~ 1 + z + x + x_second + z:x + z:x_centered_second

Misspecified:

yik ~ 1 + z + x + z:x

## Density-Ratio Estimation

InMASS density-ratio estimation uses the reconstructed x and x_second covariates. Additional squared terms are not added for this scenario.
