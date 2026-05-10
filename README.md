# metadid-sims

Simulation-based validation study for the [metadid](https://github.com/ben18785/metadid) R package. Uses the [`targets`](https://docs.ropensci.org/targets/) pipeline framework to run systematic simulation studies assessing parameter recovery, calibration, robustness to model misspecification, and edge-case behaviour.

## Quick start

### Install dependencies

metadid depends on [cmdstanr](https://mc-stan.org/cmdstanr/) and [instantiate](https://CRAN.R-project.org/package=instantiate), which compile Stan models at package install time. Install them first if you haven't already:

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()

install.packages("instantiate")
```

Then install metadid from GitHub:

```r
# install.packages("pak")
pak::pak("ben18785/metadid")
```

Install the remaining pipeline dependencies:

```r
install.packages(c("targets", "tarchetypes", "MASS", "mvtnorm", "scales", "readr"))
```

### Run the pipeline

```r
# Run the full pipeline
targets::tar_make()

# Or run individual categories
targets::tar_make(names = starts_with("A_"))  # Calibration
targets::tar_make(names = starts_with("F_"))  # Large-N bias probes
targets::tar_make(names = starts_with("B_"))  # Comparative studies
targets::tar_make(names = starts_with("C_"))  # Outlier/heavy-tailed
targets::tar_make(names = starts_with("D_"))  # Assumption violations
targets::tar_make(names = starts_with("E_"))  # Edge cases

# Render the validation report
targets::tar_make(names = "report")
```

## Configuration

Edit `_targets.R` to change the global number of replications:

```r
N_REPS <- 25L   # default; increase for tighter coverage estimates
```

All MCMC fits use `parallel_chains = 4` with 4 chains × 1000 warmup × 1000 sampling iterations.

## Scenario categories

| Category | Description | Scenarios |
|----------|------------|-----------|
| **A** | Calibration studies (coverage, bias, RMSE) | A1–A13 |
| **F** | Large-N bias probes (200 studies, narrow posteriors) | F1–F8 |
| **B** | Comparative studies (naive vs full, correlated vs independent) | B1–B4 |
| **C** | Outlier and heavy-tailed DGPs (robust vs normal) | C1–C8 |
| **D** | Assumption violations (design offsets, heterogeneous σ, misspecified ρ) | D1–D5 |
| **E** | Edge cases (extreme ρ, zero heterogeneity, unbalanced arms) | E1–E6 |

See `R/scenarios.R` for full scenario definitions.

## Outputs

### Machine-readable

- `output/aggregated_results.csv` — per-scenario summary (coverage, bias, RMSE)
- `output/replication_results.rds` — per-replication assessment data
- `_targets/` store — all intermediate results accessible via `tar_read()`

### Human-readable

- `reports/validation-report.html` — rendered Quarto report with tables and plots

## Project structure

```
metadid-sims/
├── _targets.R               # Pipeline definition
├── R/
│   ├── scenarios.R           # Scenario configurations
│   ├── simulate.R            # Simulation wrappers (metadid + bespoke DGPs)
│   ├── fit.R                 # Model fitting wrappers
│   ├── assess.R              # Assessment (coverage, bias, RMSE)
│   └── plots.R               # Summary visualisations
├── reports/
│   └── validation-report.qmd # Quarto validation report
├── output/                   # Machine-readable exports (generated)
└── README.md
```

## Adding new scenarios

1. Add a scenario definition to `SCENARIO_CONFIGS` in `R/scenarios.R`
2. Use the appropriate prefix (A, F, B, C, D, E) — the pipeline automatically discovers scenarios by prefix
3. For bespoke DGPs, add the simulation function to `R/simulate.R` and reference it via `bespoke_fn`
4. Run `tar_make()` — only new/changed targets will execute
