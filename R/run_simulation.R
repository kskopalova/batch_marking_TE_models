#===============================================================================
# run_simulation.R  -  runs a complete simulation study from the user's settings
#
#   1. checks the settings
#   2. compiles / loads the C++ files
#   3. builds the latent histories and link matrix (cached)
#   4. simulates all data sets from the true model
#   5. fits the emigration model (csv updated after every iteration; skipped if fit.TE = FALSE)
#   6. optionally fits the NoTE model to the same data (batch-marking only)
#   7. summarises the results and saves them
#===============================================================================

run.simulation <- function(model, data.type, fit.TE = TRUE, fit.noTE = FALSE,
                           n.sims, seed,
                           n.primary, n.secondary,
                           N, p, phi, alpha, beta,
                           penalty = 1, sigma = 2.5, noTE.penalty = 0,
                           root = ".") {

  #--- 1. checks ---------------------------------------------------------------
  model     <- match.arg(model, names(MODELS))
  data.type <- match.arg(data.type, c("batch", "individual"))
  if (length(n.secondary) == 1) n.secondary <- rep(n.secondary, n.primary)
  if (length(n.secondary) != n.primary)
    stop("n.secondary must have one value per primary period (n.primary = ", n.primary, ")")
  min.K <- if (model == "atbt") 4 else 3
  if (fit.TE && n.primary < min.K) stop("Model '", model, "' needs at least ", min.K, " primary periods")
  if (!fit.TE) {     # only the NoTE model is fitted: 'model' just chooses the data-generating model
    if (!model %in% c("cr", "acbc"))
      stop("With fit.TE = FALSE, model must be \"cr\" or \"acbc\" (the model the data are simulated from)")
    if (data.type == "individual") stop("fit.TE = FALSE needs batch-marking data: the NoTE model cannot be fitted to individual identification data")
    fit.noTE <- TRUE
  }
  if (fit.noTE && data.type == "individual") {
    message("NOTE: the NoTE model can only be fitted to batch-marking data; it is not fitted.")
    fit.noTE <- FALSE
  }

  S <- list(model = model, data.type = data.type, fit.TE = fit.TE, fit.noTE = fit.noTE, n.sims = n.sims, seed = seed,
            n.secondary = n.secondary, N = N, p = p, phi = phi, alpha = alpha, beta = beta,
            penalty = penalty, sigma = sigma, noTE.penalty = noTE.penalty)

  root        <- normalizePath(root)
  cpp.dir     <- file.path(root, "cpp")
  results.dir <- file.path(root, "results")
  cache.dir   <- file.path(results.dir, "cache")
  dir.create(cache.dir, showWarnings = FALSE, recursive = TRUE)

  stem <- paste0(model, "_", if (data.type == "batch") "bm" else "ind",   # e.g. atbt_bm
                 if (!fit.TE) "_noTEonly")                                # e.g. cr_bm_noTEonly
  file.emig <- file.path(results.dir, paste0(stem, "_sims.csv"))
  file.note <- file.path(results.dir, paste0(stem, if (fit.TE) "_noTE_sims.csv" else "_sims.csv"))

  #--- 2. C++ -------------------------------------------------------------------
  message("Compiling / loading C++ files...")
  load.cpp(c(if (fit.TE) get.model(model)$cpp, if (fit.noTE) "phitpt"), dir = cpp.dir)

  #--- 3. latent structure -------------------------------------------------------
  lat <- latent.structure(n.secondary, data.type, cache.dir)

  #--- 4. simulate the data ------------------------------------------------------
  message("Simulating ", n.sims, " data sets from the '", data.generating.model(model), "' model...")
  sim <- simulate.obscounts(model, lat, S)

  #--- 5. emigration model --------------------------------------------------------
  emig <- emig.summary <- NULL
  if (fit.TE) {
    emig <- run.emigration.sims(model, sim$obscounts, lat, S, file.emig)
    emig.summary <- summarise.emigration(model, emig, S, sim$true.abundance)
  }

  #--- 6. NoTE model (same data) ---------------------------------------------------
  note <- note.summary <- NULL
  if (fit.noTE) {
    str  <- noTE.structure(n.secondary, cache.dir)
    note <- run.noTE.sims(sim$obscounts, str, S, file.note, penalty = noTE.penalty)
    note.summary <- summarise.noTE(note, S, sim$true.abundance)
  }

  #--- 7. output ---------------------------------------------------------------------
  out <- list(settings = S, generated.from = sim$generated.from, obscounts = sim$obscounts, true.abundance = sim$true.abundance,
              emigration = if (fit.TE) c(emig, emig.summary),
              noTE       = if (fit.noTE) c(note, note.summary))

  ta <- sim$true.abundance                     # true abundances of every data set
  write.csv(data.frame(iteration = seq_len(n.sims),
                       setNames(as.data.frame(t(ta$avail)), paste0("Navail", seq_len(nrow(ta$avail)))),
                       setNames(as.data.frame(t(ta$alive)), paste0("Nalive", seq_len(nrow(ta$alive))))),
            file.path(results.dir, paste0(stem, "_true_abundance.csv")), row.names = FALSE)
  if (fit.TE && !is.null(emig.summary$summary))
    write.csv(emig.summary$summary, file.path(results.dir, paste0(stem, "_summary.csv")), row.names = FALSE)
  if (fit.noTE && !is.null(note.summary$summary))
    write.csv(note.summary$summary, file.path(results.dir, paste0(stem, if (fit.TE) "_noTE_summary.csv" else "_summary.csv")), row.names = FALSE)
  if (fit.noTE && !is.null(note.summary$summary.avail))   # NoTE Nalive vs true numbers available
    write.csv(note.summary$summary.avail, file.path(results.dir, paste0(stem, if (fit.TE) "_noTE_summary_avail.csv" else "_summary_avail.csv")), row.names = FALSE)

  assign(stem, out)
  save(list = stem, file = file.path(results.dir, paste0(stem, ".RData")))
  message("\nResults saved in ", results.dir, " (", stem, ".RData and csv files)")

  if (fit.TE)
    print.sim.summary(emig.summary, title = paste0("Emigration model: ", get.model(model)$label,
                      "   [data: ", if (data.type == "batch") "batch-marking" else "individual identification",
                      ", simulated from ", out$generated.from, "]"))
  if (fit.noTE) print.sim.summary(note.summary, title = paste0("NoTE model (no emigration), ",
                      if (fit.TE) "same data" else paste0("data simulated from ", out$generated.from)))

  invisible(out)
}
