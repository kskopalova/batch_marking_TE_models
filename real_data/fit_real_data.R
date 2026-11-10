
run_code <- function(){
  
  #===============================================================================
  # fit_real_data.R
  #
  # Fit the temporary-emigration models to real capture-recapture data.
  # Run this script with the working directory set to THIS folder (real_data/).
  # Needs: emigration_utils.R (in this folder) and, in ../cpp/ (shared with the simulation
  #   study): random_emigration.cpp, markovian_atbt.cpp, markovian_atbc.cpp,
  #   markovian_acbt.cpp, markovian_acbc.cpp
  #
  # Steps:  1. describe the design and enter the data
  #         2. choose models and settings
  #         3. run (compile -> prepare data -> fit -> compare -> save)
  #
  # Each model is saved as soon as it is fitted: results/<model>.RData (e.g. results/atbt.RData).
  #===============================================================================
  
  
  setwd("./")
  options(scipen=999)
  set.seed(2024)
  source("emigration_utils.R", local = TRUE)   # local = TRUE: the functions stay inside run_code() (the simulation code uses functions with the same names)
  
  #-------------------------------------------------------------------------------
  # 1. DESIGN AND DATA
  #-------------------------------------------------------------------------------
  # edit to fit your design: n.primary, n.secondary, dat
  
  n.primary   <- 4               # number of primary periods
  n.secondary <- c(2, 2, 2, 2)   # secondary occasions per primary period (one value, or one per period)
  n.secondary <- rep(n.secondary, length.out = n.primary)
  
  # observed counts as ONE numeric vector: first captures on every secondary occasion, 
  # followed by the recaptures in the order given by recapture.index(n.secondary)
  # i.e. dat <- c(...)
  dat <- "add vector of counts here"
  
  #-------------------------------------------------------------------------------
  # 2. SETTINGS
  #-------------------------------------------------------------------------------
  # which models to fit: set TRUE / FALSE
  run.model <- c(cr   = TRUE,    # completely random emigration, alpha(t)
                 atbt = TRUE,    # Markovian, alpha(t) beta(t)
                 atbc = TRUE,    # Markovian, alpha(t) beta(c)
                 acbt = TRUE,    # Markovian, alpha(c) beta(t)
                 acbc = TRUE)    # Markovian, alpha(c) beta(c)
  models   <- names(run.model)[run.model]
  penalty  <- 1                                         # 1 = penalised likelihood, 0 = none
  sigma    <- 2.5                                       # penalty strength (lambda = 1 / (2 sigma^2))
  start    <- list(alpha = 0.1, beta = 0.9, p = 0.1, phi = 0.9)   # starting values (probabilities)
  cpp.dir  <- "../cpp"                                  # folder containing the .cpp files (shared with the simulation study)
  out.dir  <- "results"                                 # all real-data outputs go here (real_data/results)
  cache    <- file.path(out.dir, "latent_structure.rds") # latent histories take a while: built once, then re-used
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  
  #-------------------------------------------------------------------------------
  # 3. RUN
  #-------------------------------------------------------------------------------
  if (is.null(dat)) stop("Enter the data vector `dat` in section 1 first.")
  
  load.models(models, dir = cpp.dir)
  
  prepared <- prepare.data(n.secondary, dat, cache = cache)
  fits     <- fit.models(models, prepared, save.dir = out.dir, penalty = penalty, sigma = sigma, start = start)
  if (length(fits) == 0) stop("None of the models could be fitted; see the messages above.")
  
  #-------------------------------------------------------------------------------
  # 4. RESULTS
  #-------------------------------------------------------------------------------
  comparison <- compare.models(fits)    # objective, AIC, BIC, N and its 95% CI
  print(comparison)
  
  fits[[1]]                             # full output of one model (here the first one); use fits$atbt, fits$cr, ...
  # fits$atbt$estimates                 # parameters with 95% CIs: N, p, phi, alpha, beta, gamma ...
  # fits$atbt$abundance                 # Navail (available for detection) and Nalive per primary period, with 95% CIs
  
  dir.create(out.dir, showWarnings = FALSE)
  saveRDS(fits, file.path(out.dir, "fits.rds"))
  write.csv(comparison, file.path(out.dir, "model_comparison.csv"), row.names = FALSE)
  for (f in fits) {
    write.csv(f$estimates, file.path(out.dir, paste0(f$model, "_parameters.csv")), row.names = FALSE)
    write.csv(f$abundance, file.path(out.dir, paste0(f$model, "_abundance.csv")),  row.names = FALSE)
  }
  
}
run_code()
