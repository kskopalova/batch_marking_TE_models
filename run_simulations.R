#===============================================================================
# run_simulations.R  -  MAIN SCRIPT
#
# Simulation study for the temporary-emigration capture-recapture models.
# Edit the settings in sections 1-4, then run the whole script.
#
# What happens:
#   * data are simulated from the 'true' model with the true values below
#       - fitted model "cr"                                -> simulated from cr
#       - fitted model atbt / atbc / acbt / acbc           -> simulated from acbc (alpha_c beta_c)
#   * the chosen emigration model is fitted to every data set
#   * optionally the NoTE model (no emigration) is fitted to the SAME data sets
#   * results are written to results/ while the simulation runs, then summarised
#===============================================================================

root <- "."          # folder that contains this script (cpp/, R/ and results/ are inside it)

#-------------------------------------------------------------------------------
# 1. WHAT TO RUN
#-------------------------------------------------------------------------------
model     <- "acbc"        # emigration model to fit: "cr", "atbt", "atbc", "acbt" or "acbc"
                           #   (if fit.TE = FALSE: the model the data are simulated from, "cr" or "acbc")
data.type <- "batch"       # "batch" = batch-marking, "individual" = individual identification
fit.TE    <- FALSE          # TRUE = fit the emigration model; FALSE = simulate only and fit the NoTE model (batch-marking only)
fit.noTE  <- TRUE         # TRUE = also fit the NoTE model (batch-marking only; ignored otherwise; always TRUE if fit.TE = FALSE)
n.sims    <- 100           # number of simulated data sets
seed      <- 2024

#-------------------------------------------------------------------------------
# 2. DESIGN
#-------------------------------------------------------------------------------
n.primary   <- 6                   # number of primary periods
n.secondary <- rep(2, n.primary)   # secondary occasions in each primary period

#-------------------------------------------------------------------------------
# 3. TRUE VALUES USED TO SIMULATE THE DATA
#-------------------------------------------------------------------------------
N     <- 5000     # superpopulation size
p     <- 0.4      # capture probability
phi   <- 0.9      # survival probability
alpha <- 0.3      # emigration probability
beta  <- 0.7      # return probability (not used when model = "cr")

#-------------------------------------------------------------------------------
# 4. FITTING OPTIONS
#-------------------------------------------------------------------------------
penalty <- 1       # emigration models: 1 = penalised likelihood, 0 = unpenalised
sigma   <- 2.5     # penalty strength (lambda = 1 / (2 * sigma^2))
noTE.penalty <- 1  # NoTE model: 1 = penalised, 0 = unpenalised

#===============================================================================
# RUN (nothing to edit below)
#===============================================================================
for (f in c("helpers.R", "true_model.R", "design.R", "emigration_models.R", "noTE.R",
            "simulate.R", "run_sims.R", "summarise.R", "run_simulation.R"))
  source(file.path(root, "R", f))

simulation <- run.simulation(model = model, data.type = data.type, fit.TE = fit.TE, fit.noTE = fit.noTE,
                             n.sims = n.sims, seed = seed,
                             n.primary = n.primary, n.secondary = n.secondary,
                             N = N, p = p, phi = phi, alpha = alpha, beta = beta,
                             penalty = penalty, sigma = sigma, noTE.penalty = noTE.penalty,
                             root = root)

# Where to find the results:
#   simulation$emigration$summary   summary table of the emigration model
#   simulation$emigration$results   raw results, one row per simulation
#   simulation$noTE$summary         the same for the NoTE model (if fitted)
#   results/<model>_<bm|ind>.RData  everything above, saved to disk
