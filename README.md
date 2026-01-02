# Integrate Meta-analysis into Specific Study (InMASS)

This repository provides an **R implementation** of **InMASS**, as introduced in the pre-print paper:  
📄 **[Integrate Meta-analysis into Specific Study (InMASS) for Estimating Conditional Average Treatment Effect](https://arxiv.org/abs/2503.21091)** 


## Repository Structure

- `01_simulation/`: Contains R scripts for conducting simulation studies.
- `02_case study/`: Contains R scripts and data for a case study.
- `propose_functions.R`: Defines the core functions used in the simulation and the case study.


## To use the core functions:

```
source("propose_functions.R")
```

## Example 
The Example of how to apply the proposed method can be found in the R scripts located in the `02_case study/` directory.




---

## Simulation Studies

### Objectives

The simulation studies are designed to assess:

- Bias and variance of the proposed estimator for the target-population treatment effect
- Type 1 error and power
- Robustness to model misspecification
- Efficiency gains from borrowing external information
- Sensitivity to sample size, number of external studies, and allocation ratios

Comparisons are made against conventional approaches such as:
- Estimation using the target trial alone
- Random-effects meta-regression–based estimators
- Partial borrowing strategies using only control-arm information

---

### Data-Generating Mechanism

Across simulations, data are generated under the following general structure:

- A **target randomized trial**, representing the population of primary interest
- Multiple **external studies**, for which only aggregated summary statistics are assumed to be available
- Covariates may differ in distribution between the target and external populations, reflecting **covariate shift**

Potential outcomes are generated under a prespecified conditional mean structure, and treatment effects are allowed to vary with covariates to induce heterogeneity in the conditional average treatment effect (CATE).

Unless otherwise stated:
- Within-study errors are independent
- Summary statistics consist of means and variances of outcomes and covariates
- External studies are treated as exchangeable conditional on covariates

---

### Estimation Procedures

For each simulated dataset, the following estimators are computed:

- The proposed **InMASS estimator**, which reconstructs second-order moments from summary statistics and applies density-ratio weighting
- Benchmark estimators based on:
  - Target trial data only
  - Meta-analytic or meta-regression approaches using external studies

Confidence intervals are constructed according to the method described in the main manuscript.

---

### Evaluation Metrics

Performance is evaluated using:

- Estimation bias
- Mean squared error
- Type 1 error and power

Results are summarized across simulation replicates and visualized using boxplots and line plots.

---

## Real-Data Analysis

### Purpose

The real-data analysis demonstrates the practical applicability of InMASS when:

- Individual-level data are available for a prespecified target population
- Only summary-level information is accessible from external studies
- Covariate distributions differ across data sources

The analysis serves as an illustrative example rather than a definitive clinical conclusion.

---

### Analysis Workflow

The real-data analysis follows these steps:

1. **Benchmark analysis**  
   A model is fitted using pooled individual-level data, when available, to obtain a reference estimate.

2. **Summary-based analysis**  
   Each external data source is treated as a data provider that shares only aggregated statistics.

3. **Application of InMASS**  
   Pseudo individual-level information is reconstructed from summary statistics and calibrated to the target population using density-ratio weighting.

4. **Comparison of results**  
   Estimates and confidence intervals from InMASS are compared with those obtained from target-only and conventional methods.

---

## Contact

For questions regarding the supplemental materials or the implementation details, please refer to the corresponding author of the manuscript.

---

## License

These materials are provided for academic and research purposes.  
Please cite the associated manuscript if you use or adapt this code.