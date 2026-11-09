#===============================================================================
# simulate.R  -  simulate observed data from the true model
#
# Data are simulated exactly as in simulate.obscounts.cr() of the original code:
# latent histories get their probabilities from the true model, a multinomial
# sample of size N is drawn over them and the link matrix turns that sample into
# the observed counts.
#
# Which model generates the data:
#     fitted model "cr"                                   -> completely random emigration (cr)
#     fitted model atbt / atbc / acbt / acbc              -> Markovian alpha(c) beta(c) (acbc)
#===============================================================================

data.generating.model <- function(model) if (model == "cr") "cr" else "acbc"

# Latent-history probabilities under the true parameter values
latent.probabilities <- function(model, lat, S) {

  secondary <- S$n.secondary
  nocca     <- sum(secondary)
  nprim     <- length(secondary)
  th        <- as.data.frame(lat$true_hists)   # the latent-probability functions expect a data.frame

  if (model == "cr") {
    a <- cr.assign(seed = S$seed, prim = nprim, occa = nocca, between.inds = secondary,
                   migrationprob = S$alpha, captureprob = S$p, survivalprob = S$phi)
    cr.probs(th = th, param_list = a$param, nOccas = nocca, sc = secondary)
  } else if (model == "acbc") {
    a <- acbc.assign(seed = S$seed, prim = nprim, occa = nocca, between.inds = secondary,
                     migrationprob = S$alpha, returnprob = S$beta,
                     captureprob = S$p, survivalprob = S$phi)
    acbc.probs(th = th, param_list = a$acbc_param, nOccas = nocca, sc = secondary)
  } else stop("No data-generating model for '", model, "'")
}

# True period-level abundances of every simulated data set.
# For each primary period the part of every latent history that belongs to that period is
# inspected: an individual is ALIVE if it is in state 1, 2 or 3 (present not captured, captured,
# away) and AVAILABLE for detection if it is in state 1 or 2. The numbers of individuals with each
# latent history in this data set are in latvec, so the counts are specific to the data set.
# Only latent histories with non-zero probability are considered.
# Returns two matrices (primary periods x simulations): alive and avail.
true.abundance.counts <- function(lat, latvec, latprobs, secondary) {

  keep <- latprobs > 0
  th   <- lat$true_hists[keep, , drop = FALSE]
  lv   <- latvec[keep, , drop = FALSE]
  pp   <- find.primperiod.indices(secondary)
  K    <- length(secondary)

  alive <- avail <- matrix(NA_real_, K, ncol(latvec))
  for (k in seq_len(K)) {
    portion  <- th[, pp[k, 1]:pp[k, 2], drop = FALSE]
    is.avail <- rowSums(portion == 1 | portion == 2) > 0
    is.alive <- rowSums(portion == 1 | portion == 2 | portion == 3) > 0
    alive[k, ] <- colSums(lv[is.alive, , drop = FALSE])
    avail[k, ] <- colSums(lv[is.avail, , drop = FALSE])
  }
  dimnames(alive) <- dimnames(avail) <- list(paste0("period", seq_len(K)), NULL)
  list(alive = alive, avail = avail)
}

# Simulate S$n.sims data sets. Returns the observed counts (one column per simulation) and the
# true abundances of every data set (determined from latvec when the data are generated).
simulate.obscounts <- function(model, lat, S) {

  gen      <- data.generating.model(model)
  latprobs <- latent.probabilities(gen, lat, S)

  set.seed(S$seed)
  latvec <- rmultinom(S$n.sims, S$N, latprobs)        # one column per simulation
  obs    <- as.matrix(lat$linkmat %*% latvec)
  dimnames(obs) <- NULL

  list(obscounts = obs, generated.from = gen, latprobs = latprobs,
       true.abundance = true.abundance.counts(lat, latvec, latprobs, S$n.secondary))
}
