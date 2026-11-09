# Simulation study: temporary-emigration capture-recapture models

```
emigration_sims/
├── run_simulations.R      MAIN SCRIPT: the only file you edit (settings) and run
├── compile_models.R       compiles the C++ files once, before the first run
├── cpp/                   all C++ files (random_emigration, markovian_*, phitpt)
├── R/                     utility scripts, sourced by the main script
│   ├── helpers.R            logit, ilogit, primary-period indices, ...
│   ├── true_model.R         true data-generating models (transition matrices, latent probabilities,
│   │                        xxx.assign / xxx.probs)  -- copied unchanged from util_emig.R
│   ├── design.R             latent histories and link matrices (batch-marking / individual id)
│   ├── simulate.R           simulate observed data from the true model
│   ├── emigration_models.R  model registry, parameter layout, compile/load, fit.model()
│   ├── noTE.R               NoTE model (phitpt): original functions + wrappers
│   ├── run_sims.R           simulation loops, csv written after every iteration
│   ├── summarise.R          summary table (true value, mean estimate, bias, coverage, CI width)
│   └── run_simulation.R     run.simulation(): the whole study from the settings
└── results/               csv, RData and summary files
```

## How to use
1. Put the C++ files in `cpp/` and run `compile_models.R` once (on the machine that will run the simulations).
2. Open `run_simulations.R`, edit sections 1-4 and run it.

## What the main script does
* **data.type** `"batch"` (batch-marking: observed data `c(m, n)`) or `"individual"` (counts of every observed capture history).
* **model** one emigration model at a time (`cr`, `atbt`, `atbc`, `acbt`, `acbc`). Data are simulated from `cr` when the fitted model is `cr`, and from `acbc` (alpha_c beta_c) for all Markovian models.
* **fit.TE** = `FALSE` skips the emigration model: data are still simulated from the true model chosen with `model` (`"cr"` or `"acbc"`) and only the NoTE model is fitted (batch-marking only). Files are then named `<model>_bm_noTEonly_...`.
* **fit.noTE** also fits the NoTE model to the *same* simulated data sets (batch-marking only; ignored for individual identification).
* Data are simulated with the true values given in the settings (latent probabilities + multinomial draw + link matrix, as in `simulate.obscounts.*`). Every fit starts from the true values.
* A csv is appended after every iteration, so a crashed run keeps what was finished.
* Results: `results/<model>_<bm|ind>_sims.csv`, `..._summary.csv` and `<model>_<bm|ind>.RData` (NoTE: `..._noTE_...`).

## Summary table
`true.value`, `mean.estimate` (not rounded), `bias`, `rel.bias.pct` (relative bias, %), `coverage` (%), `ci.width`, `n.missing` for every parameter and for the period-level abundances (`Navail1..K`, `Nalive2..`; NoTE: `Nalive1..K`).
Coverage and CI width follow `analyse.res.*` in `util_emig.R`; only simulations that converged are used.

### True period-level abundances
They are determined for every simulated data set when the data are generated (`true.abundance.counts()` in `R/simulate.R`), from the non-zero-probability latent histories and the multinomial draw `latvec` of that data set:
* **alive** in period k: individuals whose part of the latent history in period k is in state 1, 2 or 3;
* **available** for detection: state 1 or 2.

Bias = mean of (estimate - true), taken data set by data set; relative bias = mean of (estimate - true) / true x 100. For abundances `true.value` is the mean of the true values over the data sets. NoTE's `Nalive` is compared with the true numbers *alive* of the same data sets, i.e. the same reference as the emigration model. The true abundances are stored in the RData file and in `results/<model>_<bm|ind>_true_abundance.csv`.
