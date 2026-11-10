Real-data output is written here (kept apart from the simulation output in ../results):

  <model>.RData                  one fitted model (object named after the model, e.g. atbt.RData)
  <model>_parameters.csv         parameter estimates with 95% CIs
  <model>_abundance.csv          Navail and Nalive per primary period, with 95% CIs
  model_comparison.csv           objective, AIC, BIC, N and its CI for every fitted model
  fits.rds                       all fitted models
  latent_structure.rds           latent histories / link matrix for this design (slow to build, re-used)
