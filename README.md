# Latent multinomial batch-marking temporary emigration models: simulations and real data analysis

This repository contains the code for

* a **simulation study** of temporary-emigration (TE) capture-recapture models for batch-marking and individual-identification data, with an optional comparison to a model without emigration (NoTE);
* **plots** of the simulation results (ID vs batch-marking, NoTE vs TE);
* **fitting the models to real data**.

The models are written in C++ and fitted with [TMB](https://github.com/kaskr/adcomp).

---

## Contents

1. [Requirements](#requirements)
2. [Repository layout](#repository-layout)
3. [Quick start](#quick-start)
4. [Step 1: compile the models](#step-1-compile-the-models)
5. [Step 2: run simulations](#step-2-run-simulations)
6. [Step 3: make plots](#step-3-make-plots)
7. [Fitting real data](#fitting-real-data)
8. [Where everything is saved](#where-everything-is-saved)
9. [Summary tables](#summary-tables)

---

## Requirements

The code was built and run with the versions below (Windows 11 x64). Other versions may work, but these are the ones it was developed with.

| Software / package | Version | Used for |
|---|---|---|
| R | 4.5.1 | everything |
| Rtools | 4.5 (Windows only) | compiling the C++ models |
| TMB | 1.9.17, **installed from GitHub** (see below) | fitting all models |
| RcppEigen | 0.3.4.0.2 | compiling the C++ models (TMB dependency) |
| remotes | 2.5.0 | installing TMB from GitHub |
| Matrix | 1.7-3 | simulations and real data |
| tidyverse | 2.0.0 | plots |
| &nbsp;&nbsp;ggplot2 | 4.0.0 | plots |
| &nbsp;&nbsp;tidyr | 1.3.1 | plots |
| &nbsp;&nbsp;tidyselect | 1.2.1 | plots |
| scales | 1.4.0 | plots |
| patchwork | 1.3.2 | plots (combining panels) |
| rstudioapi | 0.17.1 | plots (optional: sets the working directory in RStudio) |

### Installing TMB

TMB must be installed from the GitHub fork below, **not** from CRAN:

```r
install.packages("remotes")   # if not installed
remotes::install_github("weizhangstats/adcomp/TMB")
```

TMB compiles the C++ model files, so a C++ compiler is needed:
on Windows install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version (Rtools 4.5 for R 4.5.x; `Sys.which("make")` should then return a path such as `C:\rtools45\usr\bin\make.exe`); on macOS the Xcode command-line tools; on Linux the standard build tools.

### Installing the other packages

```r
install.packages(c("Matrix", "tidyverse", "scales", "patchwork", "rstudioapi"))
```

---

## Repository layout

```
batch_marking_TE_models/
├── run_simulations.R       MAIN SCRIPT for simulations: edit the settings, then run
├── make_graphs.R           MAIN SCRIPT for plots of the simulation results
├── compile_models.R        compiles the C++ files once, before the first run
├── cpp/                    C++ model files (TMB), shared by simulations and real data
│   ├── random_emigration.cpp   completely random emigration (cr)
│   ├── markovian_atbt.cpp      Markovian emigration, alpha(t) beta(t)
│   ├── markovian_atbc.cpp      Markovian emigration, alpha(t) beta(c)
│   ├── markovian_acbt.cpp      Markovian emigration, alpha(c) beta(t)
│   ├── markovian_acbc.cpp      Markovian emigration, alpha(c) beta(c)
│   └── phitpt.cpp              NoTE model (no emigration)
├── R/                      functions, sourced by the main scripts
│   ├── helpers.R             logit, ilogit, primary-period indices, ...
│   ├── true_model.R          true data-generating models (transition matrices, latent probabilities)
│   ├── design.R              latent histories and link matrices (batch-marking / individual id)
│   ├── simulate.R            simulate observed data from the true model
│   ├── emigration_models.R   model registry, parameter layout, compile/load, fit.model()
│   ├── noTE.R                NoTE model (phitpt)
│   ├── run_sims.R            simulation loops (csv written after every iteration)
│   ├── summarise.R           summary tables (true value, mean estimate, bias, coverage, CI width)
│   ├── run_simulation.R      run.simulation(): the whole study from the settings
│   └── graph_helpers.R       loading results and plotting functions (used by make_graphs.R)
├── real_data/              SEPARATE workflow for fitting real data (not used by the simulations)
│   ├── fit_real_data.R       main script: edit the design and data, then run
│   ├── emigration_utils.R    functions for real data
│   └── results/              real-data output only
└── results/                simulation output (csv, RData)
    ├── cache/                latent histories / link matrices (created automatically, re-used)
    └── plots/                plots from make_graphs.R (created automatically)
```

`.RData`, `.csv`, `.dll` and `.o` files are excluded from git (`.gitignore`), so results and compiled models are not tracked.

---

## Quick start

1. Install the [requirements](#requirements), including TMB from GitHub.
2. Run `compile_models.R` once.
3. Edit the settings in `run_simulations.R` and run it. Results go to `results/`.
4. Edit the file names in `make_graphs.R` and run it. Plots go to `results/plots/`.

All main scripts are run with the **repository folder as the working directory** (except `real_data/fit_real_data.R`, see [below](#fitting-real-data)).

---

## Step 1: compile the models

```r
source("compile_models.R")
```

This compiles every `.cpp` file in `cpp/` and checks that it loads; a summary table is printed at the end. Compiled files (`.dll` on Windows, `.so` on Linux/macOS) are specific to the machine they were built on, so run the script on every machine that will fit models. Set `force <- FALSE` in the script to compile only the files that are missing or older than their source.
If a model has not been compiled, `run_simulations.R` and `fit_real_data.R` compile it when it is first needed.

---

## Step 2: run simulations

Open `run_simulations.R`, edit sections 1-4 and run the whole script (e.g. `source("run_simulations.R")`).

### Settings

| Section | Setting | Meaning |
|---|---|---|
| 1. What to run | `model` | emigration model to fit: `"cr"`, `"atbt"`, `"atbc"`, `"acbt"` or `"acbc"` |
| | `data.type` | `"batch"` (batch-marking) or `"individual"` (individual identification) |
| | `fit.TE` | `TRUE` = fit the emigration model; `FALSE` = only fit the NoTE model (batch-marking only) |
| | `fit.noTE` | `TRUE` = also fit the NoTE model to the **same** data sets (batch-marking only; ignored for individual identification) |
| | `n.sims`, `seed` | number of simulated data sets, random seed |
| 2. Design | `n.primary`, `n.secondary` | number of primary periods; secondary occasions in each primary period |
| 3. True values | `N`, `p`, `phi`, `alpha`, `beta` | superpopulation size, capture, survival, emigration and return probabilities (`beta` not used for `cr`) |
| 4. Fitting options | `penalty`, `sigma`, `noTE.penalty` | penalised likelihood on/off and penalty strength (lambda = 1 / (2 sigma^2)) |

### What happens

* Data are simulated from `cr` when the fitted model is `cr`, and from `acbc` (alpha_c beta_c) for all Markovian models. With `fit.TE = FALSE`, `model` is the data-generating model (`"cr"` or `"acbc"`).
* Batch-marking data are the counts `c(m, n)`; individual-identification data are the counts of every observed capture history.
* Every fit starts from the true values.
* A csv row is written after every iteration, so a crashed run keeps what was finished.
* At the end the results are summarised, printed and saved.

### Output (in `results/`)

`<model>` is the model name, `<bm|ind>` the data type (bm = batch-marking, ind = individual identification).

| File | Contents |
|---|---|
| `<model>_<bm\|ind>.RData` | **everything** (settings, simulated data, true abundances, raw results, summaries); used by the plots |
| `<model>_<bm\|ind>_sims.csv` | raw results of the emigration model, one row per simulation |
| `<model>_<bm\|ind>_summary.csv` | summary table of the emigration model |
| `<model>_bm_noTE_sims.csv`, `<model>_bm_noTE_summary.csv` | the same for the NoTE model (if fitted) |
| `<model>_bm_noTE_summary_avail.csv` | NoTE's numbers alive compared with the true numbers **available** for detection (if NoTE fitted; see [Summary tables](#note-numbers-alive-compared-with-numbers-available)) |
| `<model>_<bm\|ind>_true_abundance.csv` | true numbers alive / available of every simulated data set |
| `<model>_bm_noTEonly...` | files of a run with `fit.TE = FALSE` (NoTE model only; the NoTE files are then named `..._noTEonly_sims.csv`, `..._noTEonly_summary.csv`, `..._noTEonly_summary_avail.csv`) |
| `cache/` | latent histories and link matrices (slow to build, re-used by later runs) |

The RData file holds one object named after the file (e.g. `atbt_bm`), a list with

```
$settings         settings of the run (model, data.type, true values, ...)
$obscounts        simulated data, one column per data set
$true.abundance   true numbers $alive and $avail (rows = periods, columns = data sets)
$emigration       $results (raw), $failed, $info, $summary   (NULL if fit.TE = FALSE)
$noTE             the same for the NoTE model                (NULL if not fitted),
                  plus $summary.avail: NoTE's Nalive compared with the true numbers available
```

The same list is returned by `run.simulation()`, so after running the script it is also available as `simulation`.

**Running the same model and data type again overwrites its files.** Rename or move them first if you want to keep them.

---

## Step 3: make plots

`make_graphs.R` makes boxplots of the estimates across converged simulations, with % coverage and mean CI width printed above and below the boxes. It has two independent sections; run either one, both, or any single plot line after the loading lines of its section.

1. Copy or keep the simulation RData files in `results/`.
2. Open `make_graphs.R` and set the file names in the section you want (`list_results()` prints the files available).
3. Run the script, or the section. The working directory is set to the script's folder automatically in RStudio; otherwise run it from the repository folder.

### NoTE vs TE

Needs **one batch-marking file run with `fit.TE = TRUE` and `fit.noTE = TRUE`** (e.g. `atbt_bm.RData`). Both models were fitted to the same data sets.

```r
file_name <- "atbt_bm.RData"
results   <- load_simulation(file_name, "batch")

plot_note_vs_te("N",   results)             # superpopulation size
plot_note_vs_te("phi", results)             # survival
plot_note_vs_te_relbias("phi", results)     # survival, relative bias (%)
plot_note_vs_te_abundance("alive", results) # period-level abundance N_k, relative bias (%)
                                            #   "alive" or "available"
```

Only parameters with the same meaning in both models can be compared (`N`, `phi`, `p`), and only the periods both models estimate are shown (e.g. phi1-phi4 for `cr` and `atbt`).

### Individual identification vs batch-marking

Needs **two files for the same model**: one individual-identification run and one batch-marking run (e.g. `acbt_ind.RData` and `acbt_bm.RData`).

```r
results_id <- load_simulation("acbt_ind.RData", "individual")
results_bm <- load_simulation("acbt_bm.RData",  "batch")

plot_id_vs_bm("N", results_id, results_bm)
plot_id_vs_bm("phi", results_id, results_bm)    # also "p", "alpha", "beta", "gamma",
                                                 # "alphaprod", "betaprod", "Navail", "Nalive"
```

### Options and saving

Parameter names, the number of boxes, true values, coverage and CI widths are read from the results. Y-axis limits are set automatically; every plotting function also takes `ylim = c(lo, hi)`, `breaks = ...` and `label_y = c(cov = ..., ci = ...)` to set them by hand.

```r
save_plot(p, "atbt_note_vs_te_n", width = 6, height = 4)
```

saves the plot as png, pdf and eps (400 dpi) in **`results/plots/`** (created if needed); use `dir = ...` for another folder. Runs with `fit.TE = FALSE` (`..._noTEonly`) contain no emigration-model results and cannot be plotted.

---

## Fitting real data

Real data are fitted with the scripts in `real_data/`, which are kept apart from the simulation code (`emigration_utils.R` re-defines some function names used in `R/`, so it is sourced with `local = TRUE`). The C++ files in `cpp/` are shared with the simulations.

1. Open `real_data/fit_real_data.R`.
2. **Section 1, design and data:** set `n.primary`, `n.secondary` and the data vector `dat`. `dat` is one numeric vector of batch-marking counts: the first captures on every secondary occasion, followed by the recaptures in the order given by `recapture.index(n.secondary)`, a function in `emigration_utils.R` that returns a table of which count goes where (marking period, recapture period and recapture occasion).
3. **Section 2, settings:** choose the models to fit (`run.model`), penalty, `sigma` and starting values.
4. Run the script **with `real_data/` as the working directory**, e.g. in R

   ```r
   setwd("real_data")
   source("fit_real_data.R")
   ```

   or from a terminal: `cd real_data` then `R CMD BATCH fit_real_data.R fit_real_data.Rout`.

A model comparison table (objective, AIC, BIC, N and its 95% CI) is printed. Output goes to `real_data/results/`:

| File | Contents |
|---|---|
| `<model>.RData` | one fitted model (saved as soon as it is fitted) |
| `<model>_parameters.csv` | parameter estimates with 95% CIs |
| `<model>_abundance.csv` | numbers available (Navail) and alive (Nalive) per primary period, with 95% CIs |
| `model_comparison.csv` | objective, AIC, BIC, N and its CI for every fitted model |
| `fits.rds` | all fitted models |
| `latent_structure.rds` | latent histories / link matrix for this design (slow to build, re-used) |

Nothing in the simulation `results/` folder is touched.

---

## Where everything is saved

| What | Where |
|---|---|
| compiled models | `cpp/` |
| simulation results (RData, csv) | `results/` |
| cached latent structures (simulations) | `results/cache/` |
| plots | `results/plots/` |
| real-data results | `real_data/results/` |

---

## Summary tables

For every parameter and period-level abundance (`Navail1..K`, `Nalive2..`; NoTE: `Nalive1..K`) the summary table reports `true.value`, `mean.estimate` (not rounded), `bias`, `rel.bias.pct` (relative bias, %), `coverage` (%), `ci.width` and `n.missing`. Only simulations that converged (convergence code 0) are used; simulations with a missing estimate or interval limit are excluded from coverage.

### True period-level abundances

They are determined for every simulated data set when the data are generated (`true.abundance.counts()` in `R/simulate.R`), from the latent histories and the multinomial draw of that data set:

* **alive** in period k: individuals whose latent history in period k is in state 1, 2 or 3;
* **available** for detection: state 1 or 2.

### How the columns are calculated

Each converged simulation *i* has an estimate, interval limits Lo and Up and a true value. For parameters (N, p, phi, ...) the true value is the value used to simulate the data, the same for every data set. For period-level abundances it is the true number of **that** data set (from `$true.abundance`, matched by `iteration`), because the true numbers alive and available differ between data sets.

| Column | Calculation over the converged simulations |
|---|---|
| `true.value` | mean of the true values |
| `mean.estimate` | mean of the estimates |
| `bias` | mean of (estimate - true), taken data set by data set; equal to `mean.estimate - true.value` |
| `rel.bias.pct` | mean of (estimate - true) / true x 100, the average of the per-simulation relative errors (not `bias / true.value x 100`) |
| `coverage` | % of intervals with Lo <= true <= Up, each interval checked against the true value of its own data set; simulations with a missing estimate or limit are excluded |
| `ci.width` | mean of (Up - Lo), ignoring missing values |
| `n.missing` | number of converged simulations with a missing estimate or interval limit |

`mean.estimate`, `bias` and `rel.bias.pct` are `NA` if any converged simulation has a missing estimate (`n.missing` > 0).

NoTE's `Nalive` in `$noTE$summary` is compared with the true numbers *alive* of the same data sets, the same reference as the emigration model.

### NoTE numbers alive compared with numbers available

Without emigration, NoTE's numbers alive estimate the individuals available for detection. `$noTE$summary.avail` (and `<model>_bm_noTE_summary_avail.csv`) therefore repeats the `Nalive1..K` rows of the NoTE summary with the true numbers **available** for detection of each data set as the reference, using the same calculations as above. The estimates are the same, so `mean.estimate`, `ci.width` and `n.missing` are identical to the `Nalive` rows of `$noTE$summary`; `true.value`, `bias`, `rel.bias.pct` and `coverage` change.
