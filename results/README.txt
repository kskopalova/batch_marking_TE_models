Simulation output is written here while the simulations run:

  <model>_<bm|ind>_sims.csv          raw results, one row per simulation (updated after every iteration)
  <model>_<bm|ind>_noTE_sims.csv     same for the NoTE model (if fitted)
  <model>_<bm|ind>_summary.csv       summary table
  <model>_<bm|ind>_true_abundance.csv  true numbers alive / available of every simulated data set
  <model>_<bm|ind>_noTE_summary.csv  summary table of the NoTE model (if fitted)
  <model>_<bm|ind>.RData             everything (settings, simulated data, raw results, summaries)
  cache/                             latent histories / link matrices (re-used)

bm = batch-marking, ind = individual identification.
