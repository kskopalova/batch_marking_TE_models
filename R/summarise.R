#===============================================================================
# summarise.R  -  turn the raw simulation results into a summary table
#
# For every quantity (parameter or period-level abundance) the table reports
#   true.value     value used to simulate the data (period-level abundances: mean over the
#                  data sets of the true abundance of each data set)
#   mean.estimate  mean of the estimates over the simulations (not rounded)
#   bias           mean of (estimate - true value); the true period-level abundances differ
#                  between data sets, so the difference is taken data set by data set
#   rel.bias.pct   relative bias in %: mean of (estimate - true) / true * 100
#   coverage       % of the confidence intervals that contain the true value
#   ci.width       mean width of the confidence intervals
#   n.missing      simulations in which the estimate or an interval limit is NA/NaN
#
# Coverage and CI width are computed as in analyse.res.* of the original util_emig.R:
#   - only simulations that were fitted AND converged (convergence code 0) are used for every column
#     (failed fits have no output row; non-converged fits are removed);
#   - 'missing' = a used simulation with an NA/NaN estimate or limit (e.g. standard errors could
#     not be obtained); these are excluded from coverage, and
#     coverage = (intervals containing the truth) / (nrow(used) - n.missing) * 100;
#   - the CI width is the mean of (upper - lower) over the simulations, ignoring NA.
#===============================================================================

#-------------------------------------------------------------------------------
# TRUE VALUES
#-------------------------------------------------------------------------------

# Parameters of the fitted emigration model, in the order and with the names used in the results
true.parameters <- function(model, S) {
  K      <- length(S$n.secondary)
  nocca  <- sum(S$n.secondary)
  layout <- model.layout(model, K, nocca)
  truth  <- c(pt = S$p, phiprimt = S$phi, alphaprimt = S$alpha, alphat = S$alpha,
              betaprimt = S$beta, betat = S$beta,
              alphaprodt = S$phi * (1 - S$alpha), betaprodt = S$phi * S$beta, gamma0t = 1 / K)
  out <- c(N = S$N)
  for (b in layout$blocks)
    out <- c(out, setNames(rep(truth[[b$tmb]], b$n.rep), b$labels))
  out
}

# Names of the period-level abundances reported by the emigration model
abundance.names <- function(model, K) {
  spec    <- get.model(model)
  n.alive <- if (spec$confounded) K - 2 else K - 1
  c(paste0("Navail", seq_len(K)), paste0("Nalive", 1 + seq_len(n.alive)))
}

# True period-level abundances for the quantities named in 'nm' ("Navail<k>", "Nalive<k>"):
# a list with one vector (one value per simulated data set) per name.
# true.counts: output of true.abundance.counts() / simulate.obscounts()$true.abundance
true.abundance <- function(nm, true.counts) {
  out <- lapply(nm, function(q) {
    k <- as.integer(sub("^N(avail|alive)", "", q))
    if (startsWith(q, "Navail")) true.counts$avail[k, ] else true.counts$alive[k, ]
  })
  setNames(out, nm)
}

#-------------------------------------------------------------------------------
# SUMMARY TABLE
#-------------------------------------------------------------------------------

# res:   data.frame of raw results (columns <name>, <name>.Lo, <name>.Up for every quantity)
# truth: named list of true values; its names define the rows of the table. Each element is one
#        value (parameters) or one value per simulated data set (period-level abundances)
# type:  named vector ("parameter" / "abundance"), same names
summarise.sims <- function(res, truth, type) {

  res <- res[res$convergence == 0, , drop = FALSE]   # only simulations that converged are summarised

  one <- function(q) {
    est <- res[[q]]; lo <- res[[paste0(q, ".Lo")]]; up <- res[[paste0(q, ".Up")]]
    na.rows <- is.na(est) | is.na(lo) | is.na(up)
    tv      <- truth[[q]]
    tv      <- if (length(tv) == 1) rep(tv, nrow(res)) else tv[res$iteration]   # truth of each used data set
    n.miss  <- sum(na.rows)                       # missing = NA/NaN estimate or limit
    n.ok    <- nrow(res) - n.miss                 # converged simulations without a missing value
    data.frame(quantity = q, type = type[[q]], true.value = mean(tv),
               mean.estimate = mean(est),
               bias          = mean(est - tv),
               rel.bias.pct  = mean((est - tv) / tv) * 100,
               coverage      = if (all(is.na(tv)) || n.ok == 0) NA_real_
                               else sum(tv[!na.rows] >= lo[!na.rows] & tv[!na.rows] <= up[!na.rows]) / n.ok * 100,
               ci.width      = mean(up - lo, na.rm = TRUE),
               n.missing     = n.miss)
  }
  out <- do.call(rbind, lapply(names(truth), one))
  rownames(out) <- NULL
  out
}

# Run information: how many simulations were fitted, failed, converged
run.info <- function(res, failed, n.sims) {
  data.frame(simulations = n.sims,
             fitted      = if (is.null(res)) 0 else nrow(res),
             failed      = nrow(failed),
             not.converged  = if (is.null(res)) NA else sum(res$convergence != 0),
             used           = if (is.null(res)) 0 else sum(res$convergence == 0),
             mean.objective = if (is.null(res)) NA else mean(res$objective))
}

# Summary of an emigration-model run
summarise.emigration <- function(model, run, S, true.counts) {
  if (is.null(run$results))
    return(list(info = run.info(NULL, run$failed, S$n.sims), summary = NULL))
  K     <- length(S$n.secondary)
  par   <- as.list(true.parameters(model, S))
  ab    <- true.abundance(abundance.names(model, K), true.counts)
  truth <- c(par, ab)
  type  <- setNames(c(rep("parameter", length(par)), rep("abundance", length(ab))), names(truth))
  list(info    = run.info(run$results, run$failed, S$n.sims),
       summary = summarise.sims(run$results, truth, type))
}

# Summary of a NoTE run. Its parameters are N, p, survival phi, entry probabilities beta and
# period-level abundances (Nalive). Survival phi is compared with the phi used to simulate the
# data, exactly as for the emigration models. Without emigration the individuals alive are the
# individuals available for detection; the reference values for Nalive are the true numbers ALIVE
# of each data set, i.e. the same reference as the numbers alive of the emigration models.
summarise.noTE <- function(run, S, true.counts) {
  if (is.null(run$results))
    return(list(info = run.info(NULL, run$failed, S$n.sims), summary = NULL))
  nper  <- length(S$n.secondary)
  nocca <- sum(S$n.secondary)
  par   <- c(N = S$N,
             setNames(rep(S$p, nocca), paste0("p", seq_len(nocca))),
             setNames(rep(S$phi, nper - 1), paste0("phi", seq_len(nper - 1))),
             setNames(rep(1 / nper, nper), paste0("beta", seq_len(nper))))
  nm    <- paste0("Nalive", seq_len(nper))
  ab    <- setNames(lapply(seq_len(nper), function(k) true.counts$alive[k, ]), nm)
  truth <- c(as.list(par), ab)
  type  <- setNames(c(rep("parameter", length(par)), rep("abundance", nper)), names(truth))
  list(info    = run.info(run$results, run$failed, S$n.sims),
       summary = summarise.sims(run$results, truth, type))
}

# Readable printout (values are rounded only here, never in the stored tables)
print.sim.summary <- function(x, title = "", digits = 4) {
  op <- options(width = 160); on.exit(options(op))
  cat("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")
  cat("\nRun information\n")
  print(format(x$info, digits = digits, scientific = FALSE), row.names = FALSE)
  if (is.null(x$summary)) { cat("\nNo simulation could be fitted: see the iteration messages above.\n"); return(invisible(x)) }
  cat("\nResults (coverage and relative bias in %, NA = not available)\n")
  print(format(x$summary, digits = digits, scientific = FALSE), row.names = FALSE)
  invisible(x)
}
