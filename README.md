# Integrating Meta-analysis into a Specific Study (InMASS)

This repository contains R code for the manuscript:

**[Integrating Meta-analysis into a Specific Study (InMASS) for Estimating the Target-Population Average Treatment Effect](https://arxiv.org/abs/2503.21091)**

The manuscript has been accepted for publication in *Biometrical Journal*.

The code implements the InMASS simulation studies and the case-study analysis used in the manuscript and supplementary material.

## Repository Structure

- `01_simulation/`: Simulation-study pipeline, helper functions, scenario configuration, validation scripts, and figure-generation code.
- `01_simulation/run_final_simulations.R`: Single entry point for reproducing the full simulation workflow.
- `02_case study/`: Case-study scripts and data.
- `02_case study/rda_case-study.R`: Entry point for reproducing the case-study analysis.
- `OUTPUT_MANIFEST.csv`: Crosswalk from manuscript figures and tables to generated files, source result files, and generating scripts.
- `propose_functions.R`: Core functions used by the proposed method and related analyses.

## Requirements

The scripts require R and the standard R packages used in the simulation and case-study pipelines. Install any missing packages reported by R when running the scripts.

The simulation workflow uses parallel workers for simulation-running steps. Choose `N_WORKERS` according to the available CPU cores and memory on your machine.

## Execution Root

Run all commands in this README from the repository root: the directory that directly contains `README.md`, `OUTPUT_MANIFEST.csv`, `01_simulation/`, `02_case study/`, and `propose_functions.R`. Relative paths in the R scripts are resolved from this directory.

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

After running the simulation workflow, the nonlinear robustness MSE figure can optionally be redrawn with the target-only `3:1` method separated from the other methods:

```powershell
Rscript 01_simulation/run_plot_calibration.R
```

This post-processing script does not rerun simulations. It reads the existing nonlinear robustness results and rewrites only `robustness_nonlinear_mse.pdf` and `robustness_nonlinear_mse_3to1_only.pdf` in the nonlinear robustness figure directory.

## Recreating Figures from the Archived 10,000-Replication Results

The manuscript figures can be recreated from the 10,000-replication summary files archived on Zenodo without rerunning the simulations or extracting the entire results archive. Extract the `source_results_file` listed for the desired figure in `OUTPUT_MANIFEST.csv`, preserving its directory structure under the repository root.

The required source files are:

- `results_final_main/summary/main_results_nsim10000.csv`: Main-text Figures 2--4 and Supplementary Figures S1--S12.
- `results_final_multicov/summary/robustness_multicov_results_nsim10000.csv`: Supplementary Figures S13--S15.
- `results_final_nonlinear/summary/robustness_nonlinear_results_nsim10000.csv`: Supplementary Figures S16--S18.
- `results_final_trunc/summary/ripd_truncation_summary_all_nsim10000.csv`: Supplementary Figures S19--S20.

Run the corresponding figure-generation command from the repository root:

```powershell
Rscript 01_simulation/run_all.R --mode=make-figures --nsim=10000 --output-root=results_final_main
Rscript 01_simulation/run_all.R --mode=make-multicov-figures --nsim=10000 --output-root=results_final_multicov
Rscript 01_simulation/run_all.R --mode=make-nonlinear-figures --nsim=10000 --output-root=results_final_nonlinear
Rscript 01_simulation/run_all.R --mode=make-ripd-truncation-figures --nsim=10000 --output-root=results_final_trunc
```

The `--output-root` argument identifies both the directory containing the input summary file and the directory in which the recreated figures are written. The `--nsim=10000` argument must be specified explicitly because the default for the figure-generation modes is 10. Each command recreates all figures in the corresponding simulation family rather than a single PDF in isolation.

To recreate the two nonlinear MSE components used together in Supplementary Figure S17, first run the nonlinear figure command above, set `OUTPUT_ROOT <- "results_final_nonlinear"` near the top of `01_simulation/run_plot_calibration.R`, and then run:

```powershell
Rscript 01_simulation/run_plot_calibration.R
```

## Reproducing the Case Study

The case-study analysis can be reproduced by running:

```powershell
Rscript "02_case study/rda_case-study.R"
```

The path contains a space, so quotes are recommended when running the command from a shell.

The case-study script writes descriptively named diagnostic plots and the CSV files used to reproduce Main-text Tables 1--3 to `02_case study/`.

## Output-to-Manuscript Crosswalk

`OUTPUT_MANIFEST.csv` maps Main-text Figures 2--4 and Tables 1--3 and Supplementary Figures S1--S20 to their generated output files, source result files, and generating scripts. For quick test runs with `NSIM = 10`, replace the `nsim10000` portion of a source result filename with the corresponding test-run value.

## Using Core Functions

To load the core functions directly in an R session:

```r
source("propose_functions.R")
```

## Archived Results and Version Citation

The complete 10,000-replication simulation results and logs will be archived on Zenodo:

- Zenodo DOI: `doi:xxx` (placeholder)

For exact reproducibility, cite the Zenodo results record together with a fixed GitHub Release and the full commit hash used for that Release. Do not cite only the moving `main` branch. The final citation information will be provided in the following form:

- GitHub Release: `https://github.com/keisuke-hanada/InMASS/releases/tag/<release-tag>`
- Git commit: `https://github.com/keisuke-hanada/InMASS/commit/<full-commit-hash>`

## Contact

For questions regarding the supplemental materials or the implementation details, please refer to the corresponding author of the manuscript.

## License

These materials are provided for academic and research purposes. Please cite the associated manuscript if you use or adapt this code.
