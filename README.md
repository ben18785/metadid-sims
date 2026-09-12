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

**The Quarto CLI must be on `PATH`**, even to run a single category.
`tarchetypes::tar_quarto()` probes the report at pipeline-CONSTRUCTION time, so
without it `tar_make()` fails for every target with `Quarto CLI not found` — not
just for the report.

```r
# Run the full pipeline
targets::tar_make()

# Or run individual categories
targets::tar_make(names = starts_with("A_"))  # Calibration
targets::tar_make(names = starts_with("X_"))  # Randomisation / baseline imbalance

# Render the validation report
targets::tar_make(names = "report")
```

The pipeline fits against the **installed** `metadid`, while `_targets.R` tracks
`../metadid/R` for invalidation. Those can disagree: editing metadid invalidates
every target, which then re-runs against an unchanged installed build and
attributes the results to the new code. Reinstall before running:

```sh
R CMD INSTALL ../metadid
```

### Run the unit tests

Fast, deterministic checks of the simulation scaffolding. They fit nothing and
finish in seconds, so they are safe to run on every change:

```r
testthat::test_dir("tests/testthat")
```

See [Tests](#tests) for what they cover and why.

## Configuration

Edit `_targets.R` to change the global number of replications:

```r
N_REPS <- 25L   # default; increase for tighter coverage estimates
```

All MCMC fits use `parallel_chains = 4` with 4 chains × 1000 warmup × 1000 sampling iterations.

## Scenario categories

256 scenarios across 24 categories, expanding to 460 fitted models once
per-scenario comparison arms are counted.

**Validation categories** — is the model calibrated, and where does it break?

| Category | Description | n |
|----------|-------------|--:|
| **A** | Calibration (coverage, bias, RMSE), plus sign and direction variants | 18 |
| **B** | Comparative studies (naive vs full, correlated vs independent) | 6 |
| **C** | Outlier and heavy-tailed DGPs (robust vs normal) | 8 |
| **D** | Assumption violations (design offsets, heterogeneous σ, misspecified ρ) | 5 |
| **E** | Edge cases (extreme ρ, zero heterogeneity, unbalanced arms) | 6 |
| **F** | Large-N bias probes (200 studies, narrow posteriors) | 8 |
| **G** | Bias source investigation | 9 |
| **H** | Time-trend distributional misspecification | 4 |
| **I** | Multiplicative covariates | 9 |
| **X** | **Randomisation and baseline imbalance** | 27 |

**Figure sweeps** — parameter sweeps feeding the paper's panels. These run at a
reduced replication count (`N_REPS_FIG`) and are built by `figure.yaml` rather
than the daily pipeline.

| Category | Description | n |
|----------|-------------|--:|
| **J** | Trend-mean sweep (panel B) | 9 |
| **K** | Composition sweep, small DiD core (panel B) | 17 |
| **L** | Outlier-fraction sweep (panel E) | 5 |
| **M** | Trend-variability sweep (panel C) | 6 |
| **N** | Exchangeability-gap sweep (panel F) | 7 |
| **O** | Trend-crossover sweep (panel G) | 9 |
| **P** | Unrandomised post-only imbalance sweep | 7 |
| **Q** | Trend-plane grid (panels K and L) | 25 |
| **R** | How many DiD studies buy back the premium (panel J) | 18 |
| **S** | Incomplete designs under high effect heterogeneity | 9 |
| **T** | M sweep at a PP-heavy composition (panel F) | 6 |
| **U** | What each incomplete design is worth (panel E) | 24 |
| **V** | K composition sweep at the S core size (panel B) | 9 |
| **W** | Does the naive weight share depend on trend heterogeneity? | 5 |

Category **X** asks whether baseline imbalance learned from one kind of study
should transport to another — the question behind metadid's two-population
`gamma` model. It covers transport bias, direction vs magnitude, the cost of a
wrong randomisation label, the hard-zero trap, `kappa` identification, cluster
randomisation, design-offset aliasing, the normalisation divisor, individual
data, and the pre-post scale mismatch. See the category X section of the
validation report.

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
│   ├── scenarios.R           # Scenario configurations and expectations registry
│   ├── simulate.R            # Simulation wrappers (metadid + bespoke DGPs)
│   ├── fit.R                 # Model fitting wrappers
│   ├── fit_g.R               # Category G fitting variant
│   ├── fit_naive.R           # Naive comparator arm
│   ├── fit_robust.R          # Robust-heterogeneity comparator arm
│   ├── assess.R              # Assessment (coverage, bias, RMSE)
│   ├── illustration.R        # Worked illustration for the paper
│   └── plots.R               # Summary visualisations
├── tests/testthat/           # Fast unit tests (fit nothing, run in seconds)
├── reports/
│   └── validation-report.qmd # Quarto validation report
├── docs/                     # Derivations supporting the figure panels
├── output/                   # Machine-readable exports (generated)
├── test_smoke.R              # One replication per category, end to end (slow)
├── validate_estimand.R       # Focused estimand checks (slow)
├── validate_g_subset.R       # Focused category G checks (slow)
└── README.md
```

## Tests

Two layers, with different purposes.

**`tests/testthat/` — fast, deterministic, fits nothing.** Runs in seconds on
every push via `unit-tests.yaml`. In a simulation study the expensive part is
the MCMC, but the dangerous part is the deterministic scaffolding around it: a
wrong truth value, a mis-mapped parameter, or a DGP that does not generate what
the scenario asked for produces a plausible wrong number rather than an error.

| File | Covers |
|------|--------|
| `test-dgp-moments.R` | simulated data carries the θ, β, baselines, σ and ρ it was configured with — checked from raw means, independently of `build_true_params()` |
| `test-true-params.R` | the delta-method truths every bias and coverage number is scored against |
| `test-simulate-imbalance.R` | randomisation DGP invariants — labels, mislabelling, clustering, γ populations |
| `test-assess.R` | scale selection, pinned parameters dropped rather than credited |
| `test-scoring-completeness.R` | every reported parameter is scored or explicitly excused |
| `test-scenario-config.R` | simulators resolve, arms are labelled, and every fit-config key is a real argument of the metadid function it targets |

That last check is what catches metadid's API drifting out from under the sims —
in milliseconds, rather than inside a fit.

**`test_smoke.R` — one replication per category, end to end (slow).** Runs the
real config → simulate → fit → assess path at minimal MCMC settings. Run it by
hand after changing the pipeline's plumbing.

## Adding new scenarios

1. Add a scenario definition to `SCENARIO_CONFIGS` in `R/scenarios.R`
2. Use an existing category prefix — the pipeline discovers scenarios by prefix
   via `scenario_ids()`, and a new prefix also needs a target block in
   `_targets.R` and an entry in `all_agg` / `all_rep`
3. For bespoke DGPs, add the simulation function to `R/simulate.R` and reference it via `bespoke_fn`
4. Run `tar_make()` — only new/changed targets will execute

## Acknowledgements

Parts of this simulation study were developed with assistance from
[Claude Code](https://claude.com/claude-code) (Anthropic), used as a tool under
the direction and review of the authors.
