#===============================================================================
# run_sims.R  -  the simulation loops
#
# Every iteration is appended to a csv file as soon as it has been fitted, so if a
# run crashes the results of the finished iterations are kept.
#===============================================================================

# Append one row to a csv (writes the header with the first row)
append.csv <- function(row, file) {
  first <- !file.exists(file)
  write.table(row, file, sep = ",", row.names = FALSE, col.names = first, append = !first)
}

# Generic loop: fit.once(i) must return a one-row data.frame (or fail with an error)
sim.loop <- function(n.sims, label, csv, fit.once) {
  if (file.exists(csv)) file.remove(csv)
  rows <- vector("list", n.sims)
  failed <- data.frame(iteration = integer(0), error = character(0))
  for (i in seq_len(n.sims)) {
    message("\n[", label, "]  iteration ", i, " of ", n.sims, "   (", format(Sys.time(), "%H:%M:%S"), ")")
    row <- tryCatch(suppressWarnings(fit.once(i)), error = function(e) {
      message("  iteration ", i, " failed: ", conditionMessage(e))
      failed[nrow(failed) + 1, ] <<- list(i, conditionMessage(e))
      NULL
    })
    if (!is.null(row)) {
      rows[[i]] <- row
      append.csv(row, csv)
    }
  }
  res <- do.call(rbind, rows)
  rownames(res) <- NULL
  list(results = res, failed = failed)
}

#-------------------------------------------------------------------------------
# Emigration model
#-------------------------------------------------------------------------------

# Flatten a fitted emigration model into one csv row:
# iteration, objective, convergence | estimates | lower limits (.Lo) | upper limits (.Up)
# Abundance: Navail1..NavailK, then Nalive2.. (Nalive1 = Navail1; the last Nalive is missing
# for models with a confounded final period).
emigration.row <- function(fit, iteration) {
  e   <- fit$estimates
  ab  <- fit$abundance
  al  <- ab[ab$period >= 2 & !is.na(ab$Nalive), ]
  nm  <- c(e$parameter, paste0("Navail", ab$period), paste0("Nalive", al$period))
  est <- c(e$estimate, ab$Navail,    al$Nalive)
  lo  <- c(e$lower,    ab$Navail.lo, al$Nalive.lo)
  up  <- c(e$upper,    ab$Navail.up, al$Nalive.up)
  out <- c(iteration = iteration, objective = fit$nllike, convergence = fit$convergence,
           setNames(est, nm), setNames(lo, paste0(nm, ".Lo")), setNames(up, paste0(nm, ".Up")))
  as.data.frame(as.list(out), check.names = FALSE)
}

# Minimum population size = number of distinct individuals captured.
#   batch-marking:  the first nocca counts (m) are the first captures on each occasion
#   individual id:  every observed history is a distinct individual, so the individuals first
#                   captured on each occasion add up to the sum of all observed counts
min.population <- function(y, nocca, data.type)
  if (data.type == "batch") sum(y[seq_len(nocca)]) else sum(y)

# Fit one emigration model to every simulated data set. The starting values are the true values.
run.emigration.sims <- function(model, obscounts, lat, S, csv) {

  K     <- length(S$n.secondary)
  nocca <- sum(S$n.secondary)
  start <- list(alpha = S$alpha, beta = S$beta, p = S$p, phi = S$phi)

  fit.once <- function(i) {
    y   <- obscounts[, i]
    dat <- list(nprim = K, nocca = nocca, ppInd = find.primperiod.indices(S$n.secondary),
                true_hists = lat$true_hists, linkmat = lat$linkmat,
                obscounts = y, Nmin = min.population(y, nocca, S$data.type))
    fit <- fit.model(model, dat, penalty = S$penalty, sigma = S$sigma, start = start,
                     Nt.start = log(S$N - dat$Nmin))
    emigration.row(fit, i)
  }
  sim.loop(S$n.sims, paste0("emigration model ", model), csv, fit.once)
}

#-------------------------------------------------------------------------------
# NoTE model
#-------------------------------------------------------------------------------
run.noTE.sims <- function(obscounts, str, S, csv, penalty = 0) {
  fit.once <- function(i) fit.noTE.once(obscounts[, i], i, str, S, penalty)
  sim.loop(S$n.sims, "NoTE model", csv, fit.once)
}
