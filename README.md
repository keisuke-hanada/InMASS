# Integrating Meta-analysis into a Specific Study (InMASS)

This repository contains R code for the manuscript:

**[Integrating Meta-analysis into a Specific Study (InMASS) for Estimating the Target-Population Average Treatment Effect](https://arxiv.org/abs/2503.21091)**

The code implements the InMASS simulation studies and the case-study analysis used in the manuscript and supplementary material.

## Repository Structure

- `01_simulation/`: Simulation-study pipeline, helper functions, scenario configuration, validation scripts, and figure-generation code.
- `01_simulation/run_final_simulations.R`: Single entry point for reproducing the full simulation workflow.
- `02_case study/`: Case-study scripts and data.
- `02_case study/rda_case-study.R`: Entry point for reproducing the case-study analysis.
- `propose_functions.R`: Core functions used by the proposed method and related analyses.

## Requirements

The scripts require R and the standard R packages used in the simulation and case-study pipelines. Install any missing packages reported by R when running the scripts.

The simulation workflow uses parallel workers for simulation-running steps. Choose `N_WORKERS` according to the available CPU cores and memory on your machine.

## Reproducing the Simulation Studies

The full simulation workflow is controlled by:

```r
01_simulation/run_final_simulations.R
```

At the top of that script, edit the user-configurable constants:

```r
NSIM <- 10
N_WORKERS <- 20
```

Use `NSIM <- 10` for a quick test run. For the manuscript-level simulation run, set:

```r
NSIM <- 10000
```

and adjust `N_WORKERS` as appropriate for your machine.

Run the complete workflow from the repository root with:

```powershell
Rscript 01_simulation/run_final_simulations.R
```

The script runs all simulation components sequentially:

- main one-covariate simulation scenarios;
- multi-covariate robustness scenarios;
- nonlinear robustness scenarios;
- RIPD residual-variance truncation diagnostics;
- internal validation steps;
- figure generation.

The script writes results, summaries, figures, logs, and a run manifest to the output roots configured near the top of `run_final_simulations.R`, such as:

- `results_final_main/`
- `results_final_multicov/`
- `results_final_nonlinear/`
- `results_final_trunc/`
- `logs_final/final_simulation_manifest.csv`

## Reproducing the Case Study

The case-study analysis can be reproduced by running:

```powershell
Rscript "02_case study/rda_case-study.R"
```

The path contains a space, so quotes are recommended when running the command from a shell.

## Using Core Functions

To load the core functions directly in an R session:

```r
source("propose_functions.R")
```

## Contact

For questions regarding the supplemental materials or the implementation details, please refer to the corresponding author of the manuscript.

## License

These materials are provided for academic and research purposes. Please cite the associated manuscript if you use or adapt this code.
