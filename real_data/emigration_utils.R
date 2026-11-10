#===============================================================================
# emigration_utils.R
#
# Helper functions for fitting the temporary-emigration models to real
# capture-recapture data (use together with fit_real_data.R).
#
# Models (the name is what you pass to the fitting functions):
#   "cr"   completely random emigration, time-dependent alpha
#   "atbt" Markovian emigration, time-dependent alpha, time-dependent beta
#   "atbc" Markovian emigration, time-dependent alpha, constant beta
#   "acbt" Markovian emigration, constant alpha,       time-dependent beta
#   "acbc" Markovian emigration, constant alpha,       constant beta
#
# Main functions:
#   recapture.index()  order in which the recapture counts n must be supplied
#   prepare.data()     data -> latent histories, link matrix, observed counts
#   load.models()      compile / load the C++ files
#   fit.model()        fit ONE model
#   fit.models()       fit several models
#   compare.models()   one-row-per-model summary table
#===============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(TMB)
})

#-------------------------------------------------------------------------------
# 1. MODEL REGISTRY
#    The only place where the models differ. "confounded" models are those where
#    the last survival/emigration (and return) parameters cannot be separated,
#    so only the products phi*(1-alpha) (and phi*beta) are estimated.
#-------------------------------------------------------------------------------
MODELS <- list(
  cr   = list(cpp = "random_emigration", alpha = "time",     beta = "none",     confounded = TRUE,
              label = "Completely random emigration, alpha(t)"),
  atbt = list(cpp = "markovian_atbt",    alpha = "time",     beta = "time",     confounded = TRUE,
              label = "Markovian emigration, alpha(t) beta(t)"),
  atbc = list(cpp = "markovian_atbc",    alpha = "time",     beta = "constant", confounded = FALSE,
              label = "Markovian emigration, alpha(t) beta(c)"),
  acbt = list(cpp = "markovian_acbt",    alpha = "constant", beta = "time",     confounded = FALSE,
              label = "Markovian emigration, alpha(c) beta(t)"),
  acbc = list(cpp = "markovian_acbc",    alpha = "constant", beta = "constant", confounded = FALSE,
              label = "Markovian emigration, alpha(c) beta(c)")
)

get.model <- function(model) {
  if (!model %in% names(MODELS))
    stop("Unknown model '", model, "'. Choose from: ", paste(names(MODELS), collapse = ", "))
  MODELS[[model]]
}

#-------------------------------------------------------------------------------
# 2. SMALL HELPERS
#-------------------------------------------------------------------------------
logit  <- function(x) log(x / (1 - x))
ilogit <- function(x) 1 / (1 + exp(-x))

# Entry probabilities -> transformed 'conditional entry' parameters used by the C++ code
trans.beta.func <- function(beta0) {
  out <- numeric(length(beta0))
  cs  <- 0
  for (i in seq_along(beta0)) {
    out[i] <- logit(beta0[i] / (1 - cs))
    cs <- cs + beta0[i]
  }
  out
}

# First and last secondary occasion of every primary period
find.primperiod.indices <- function(secondary) {
  last  <- cumsum(secondary)
  first <- last - secondary + 1
  ind <- cbind(first, last)
  dimnames(ind) <- list(paste("Period", seq_along(secondary)),
                        c("Starting occasion", "Ending occasion"))
  ind
}

#-------------------------------------------------------------------------------
# 3. DATA y = (m, n)
#    y = c(m, n) is a single vector of counts:
#    m: first captures on each secondary occasion        (length = total occasions)
#    n: recaptures of previously marked individuals       (see recapture.index)
#-------------------------------------------------------------------------------

# Table describing every element of n, in the order they must be supplied:
# for each primary period in which animals were marked (MarkPrimary), the recapture
# counts on every later occasion of that period (excluding its first occasion) and
# on all occasions of all later primary periods.
recapture.index <- function(secondary) {
  K   <- length(secondary)
  occ <- data.frame(Primary   = rep(seq_len(K), secondary),
                    Secondary = sequence(secondary))
  do.call(rbind, lapply(seq_len(K), function(p1) {
    later <- occ[occ$Primary > p1 | (occ$Primary == p1 & occ$Secondary != 1), ]
    data.frame(MarkPrimary = p1, RecapPrimary = later$Primary, RecapSecondary = later$Secondary)
  }))
}

#-------------------------------------------------------------------------------
# 4. LATENT HISTORIES AND LINK MATRIX (depend on the design only, not on the data)
#    Latent states: 0 not yet entered, 1 present/not captured, 2 first capture,
#                   3 temporarily away, 4 dead.
#-------------------------------------------------------------------------------
generate.latent.histories <- function(sc) {

  nprim <- length(sc)
  nocca <- sum(sc)
  ppind <- find.primperiod.indices(sc)

  # Step 1: skeleton histories (0 = not entered, 6 = alive, 4 = dead)
  i.ind <- c(0, cumsum(sc)[-length(sc)])
  j.ind <- rev(c(0, cumsum(rev(sc)[-length(sc)])))
  skeleton <- NULL
  for (i in seq_along(i.ind)) {
    for (j in seq(i, length(j.ind))) {
      mid <- nocca - i.ind[i] - j.ind[j]
      grid <- do.call(expand.grid, as.list(c(rep(0, i.ind[i]), rep(6, mid), rep(4, j.ind[j]))))
      skeleton <- rbind(skeleton, as.matrix(grid))
    }
  }
  colnames(skeleton) <- NULL

  first.of.period <- function(h) apply(ppind, 1, function(ind) h[ind[1]:ind[2]][1])

  replace.six <- function(vec, comb) { vec[vec == 6] <- comb; vec }

  # Step 2: decide, period by period, whether an alive animal is in the study area (6) or away (3)
  expand.away <- function(row) {
    th <- first.of.period(row)
    zeros <- which(th == 0)
    fours <- which(th == 4)
    last.zero  <- if (length(zeros)) tail(zeros, 1) else 0
    first.four <- if (length(fours)) head(fours, 1) else nprim + 1
    middle <- (last.zero + 1):(first.four - 1)
    combs  <- expand.grid(rep(list(c(3, 6)), length(middle)))
    res    <- t(apply(combs, 1, replace.six, vec = th))
    # repeat each period's state on all its secondary occasions
    res[, rep(seq_len(nprim), sc), drop = FALSE]
  }
  comb <- do.call(rbind, lapply(seq_len(nrow(skeleton)), function(i) expand.away(skeleton[i, ])))
  comb <- comb[comb[, 1] != 3, , drop = FALSE]
  comb <- comb[rowSums(comb == 6) > 0, , drop = FALSE]

  # Step 3: in every 'present' (6) occasion the animal is either not captured (1) or captured (2)
  expand.capture <- function(row) {
    combs <- expand.grid(rep(list(c(1, 2)), sum(row == 6)))
    t(apply(combs, 1, replace.six, vec = row))
  }
  th <- do.call(rbind, lapply(seq_len(nrow(comb)), function(i) expand.capture(comb[i, ])))
  dimnames(th) <- NULL
  th
}

# Link matrix: rows = observed counts (c(m, n) in the order of recapture.index),
# columns = latent histories. Returned as a sparse matrix.
find.linkmat <- function(th, sc) {
  ppInd <- find.primperiod.indices(sc)
  nhist <- nrow(th)
  n     <- ncol(th)
  nper  <- nrow(ppInd)

  # recaptures = every capture after the first one (state 2 -> 5)
  th5 <- th
  for (i in seq_len(nhist)) {
    cap <- which(th5[i, ] == 2)
    if (length(cap) > 1) th5[i, cap[-1]] <- 5
  }

  # A: marked individuals, one row per occasion
  Amat <- matrix(0, nrow = n, ncol = nhist)
  for (i in seq_len(nhist)) Amat[which(th5[i, ] == 2), i] <- 1

  # B: recaptured individuals, one block of n rows per primary period of marking
  Bmat <- matrix(0, nrow = nper * n, ncol = nhist)
  for (i in seq_len(nhist)) {
    if (any(th5[i, ] == 5)) {
      markocca <- which(th5[i, ] == 2)
      markpp   <- which((ppInd[, 1] - markocca) * (ppInd[, 2] - markocca) <= 0)
      for (j in ((markpp - 1) * n + markocca + 1):(markpp * n)) {
        if (th5[i, j - n * (markpp - 1)] == 5) Bmat[j, i] <- 1
      }
    }
  }
  Bmat <- Bmat[rowSums(Bmat) > 0, , drop = FALSE]   # rows that can never be observed

  Matrix(rbind(Amat, Bmat), sparse = TRUE)
}

#-------------------------------------------------------------------------------
# 5. PREPARE THE DATA
#    cache: optional .rds file; the histories/link matrix are slow to build and only
#           depend on 'secondary', so they are stored and re-used on later runs.
#-------------------------------------------------------------------------------
prepare.data <- function(secondary, y, cache = NULL) {

  K     <- length(secondary)
  nocca <- sum(secondary)

  lat <- NULL
  if (!is.null(cache) && file.exists(cache)) {
    lat <- readRDS(cache)
    if (!identical(as.numeric(lat$secondary), as.numeric(secondary))) lat <- NULL
  }
  if (is.null(lat)) {
    message("Generating latent histories and link matrix (slow, done only once if 'cache' is set)...")
    th  <- generate.latent.histories(secondary)
    lat <- list(secondary = secondary, true_hists = th, linkmat = find.linkmat(th, secondary))
    if (!is.null(cache)) saveRDS(lat, cache)
  }

  structure(list(secondary = secondary, nprim = K, nocca = nocca,
                 ppInd = find.primperiod.indices(secondary),
                 true_hists = lat$true_hists, linkmat = lat$linkmat,
                 obscounts = y, Nmin = sum(y[seq_len(nocca)])),
            class = "emig_data")
}

# Observed counts equal to zero are removed together with the latent histories that generate them
drop.zero.counts <- function(dat) {
  zero <- which(dat$obscounts == 0)
  if (length(zero) == 0) return(dat)
  bad <- colSums(dat$linkmat[zero, , drop = FALSE]) > 0
  dat$linkmat    <- dat$linkmat[-zero, !bad, drop = FALSE]
  dat$true_hists <- dat$true_hists[!bad, , drop = FALSE]
  dat$obscounts  <- dat$obscounts[-zero]
  dat
}

#-------------------------------------------------------------------------------
# 6. MODEL LAYOUT: which parameters each model has, how many, and their starting values.
#    Blocks are listed in the order in which the C++ code reports them (ADREPORT),
#    followed by N_avail (K values) and N_alive (K-1 values, K-2 if confounded).
#-------------------------------------------------------------------------------
default.start <- function() list(alpha = 0.1, beta = 0.9, p = 0.1, phi = 0.9)

model.layout <- function(model, K, nocca, start = default.start()) {
  spec <- get.model(model)
  if (K < (if (model == "atbt") 4 else 3))
    stop("Model '", model, "' needs at least ", if (model == "atbt") 4 else 3, " primary periods")

  nT <- if (spec$confounded) K - 2 else K - 1   # number of free phi (and time-dependent alpha)
  nB <- if (spec$confounded) K - 3 else K - 2   # number of free time-dependent beta

  # one block = one TMB parameter. start is given on the probability scale.
  block <- function(label, tmb, n, prob, first = 1, constant = FALSE)
    list(tmb = tmb, n = n, n.rep = n, start = rep(logit(prob), length.out = n),
         labels = if (constant) label else paste0(label, seq(first, length.out = n)),
         link = "logit")

  blocks <- list(block("p", "pt", nocca, start$p),
                 block("phi", "phiprimt", nT, start$phi))
  blocks <- c(blocks, list(
    if (spec$alpha == "time") block("alpha", "alphaprimt", nT, start$alpha)
    else                      block("alpha", "alphat", 1, start$alpha, constant = TRUE)))
  if (spec$beta == "time")     blocks <- c(blocks, list(block("beta", "betaprimt", nB, start$beta, first = 2)))
  if (spec$beta == "constant") blocks <- c(blocks, list(block("beta", "betat", 1, start$beta, constant = TRUE)))
  if (spec$confounded) {
    blocks <- c(blocks, list(block("alphaprod", "alphaprodt", 1, start$phi * (1 - start$alpha), constant = TRUE)))
    if (spec$beta != "none")
      blocks <- c(blocks, list(block("betaprod", "betaprodt", 1, start$phi * start$beta, constant = TRUE)))
  }
  # entry probabilities: K-1 free parameters, K probabilities reported
  blocks <- c(blocks, list(list(tmb = "gamma0t", n = K - 1, n.rep = K,
                                start = trans.beta.func(rep(1 / K, K - 1)),
                                labels = paste0("gamma", seq_len(K)), link = "logit")))

  list(spec = spec, blocks = blocks,
       n.avail = K, n.alive = if (spec$confounded) K - 2 else K - 1)
}

# List of everything the C++ code declares as DATA
tmb.data <- function(layout, dat, penalty, sigma) {
  spec <- layout$spec
  K    <- dat$nprim
  pp   <- dat$ppInd[, 2]
  lam  <- 1 / (2 * sigma^2)

  out <- list(obs_counts = dat$obscounts,
              linkmat    = as(dat$linkmat, "TsparseMatrix"),
              true_hists = as.matrix(dat$true_hists),
              Nmin  = dat$Nmin, nocca = dat$nocca, nhist = nrow(dat$true_hists), nper = K,
              primindend   = pp[-K],
              primindstart = dat$ppInd[, 1],
              phi_lambda = rep(lam, K - 1), p_lambda = rep(lam, dat$nocca),
              gamma0_lambda = rep(lam, K - 1),
              a_lambda = rep(lam, if (spec$alpha == "time") K - 1 else 1),
              penalty = penalty)

  if (spec$confounded) {
    out$primindend_par   <- pp[-c(K - 1, K)]
    out$alphaprod_lambda <- lam
  }
  if (spec$beta != "none") {
    out$subprimindend <- if (spec$confounded) pp[-c(K - 1, K)][-1] else pp[-K][-1]
    out$b_lambda      <- rep(lam, if (spec$beta == "time") K - 2 else 1)
    if (spec$confounded) out$betaprod_lambda <- lam
  }
  out
}

#-------------------------------------------------------------------------------
# 7. COMPILE / LOAD THE C++ FILES
#-------------------------------------------------------------------------------
load.models <- function(models = names(MODELS), dir = ".", flags = "-O1 -g") {
  old <- setwd(dir); on.exit(setwd(old))
  for (mod in models) {
    cpp <- get.model(mod)$cpp
    src <- paste0(cpp, ".cpp")
    if (!file.exists(dynlib(cpp)) || file.mtime(dynlib(cpp)) < file.mtime(src))
      compile(src, flags, DLLFLAGS = "")
    if (!cpp %in% names(getLoadedDLLs())) dyn.load(dynlib(cpp))
  }
  invisible(models)
}

#-------------------------------------------------------------------------------
# 8. FIT ONE MODEL
#-------------------------------------------------------------------------------
# dat       output of prepare.data()
# penalty   1 = penalised likelihood (ridge penalty on the logit parameters), 0 = none
# sigma     penalty strength: lambda = 1 / (2 sigma^2)
# start     starting values (probabilities) for p, phi, alpha, beta
# Nt.start  starting value for log(N - Nmin)
# abundance.ci  "wald" (symmetric, count scale) or "log" (log-normal) intervals for Navail / Nalive
fit.model <- function(model, dat, penalty = 1, sigma = 2.5, start = default.start(),
                      Nt.start = 1, abundance.ci = c("wald", "log"), silent = TRUE) {

  abundance.ci <- match.arg(abundance.ci)
  z       <- qnorm(0.975)
  K       <- dat$nprim
  layout  <- model.layout(model, K, dat$nocca, start)
  dat     <- drop.zero.counts(dat)

  # --- TMB object and optimisation
  pars <- c(list(Nt = Nt.start),
            setNames(lapply(layout$blocks, `[[`, "start"), sapply(layout$blocks, `[[`, "tmb")),
            list(svec = rep(1, length(dat$obscounts))))
  obj <- MakeADFun(tmb.data(layout, dat, penalty, sigma), pars, random = "svec",
                   DLL = layout$spec$cpp, silent = silent, saddlepoint = TRUE, hessian = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 3000, eval.max = 1000))
  if (opt$convergence != 0) warning("Model '", model, "': optimiser did not converge (", opt$message, ")")

  sdr <- sdreport(obj)
  raw <- summary(sdr, select = "report")
  if (!isTRUE(sdr$pdHess)) warning("Model '", model, "': Hessian is not positive definite; standard errors are unreliable")

  # --- split the reported quantities: Nt | parameter blocks | Navail | Nalive
  n.rep.par <- sum(sapply(layout$blocks, `[[`, "n.rep"))
  n.nk      <- layout$n.avail + layout$n.alive
  if (nrow(raw) != 1 + n.rep.par + n.nk)
    stop("Model '", model, "': C++ reports ", nrow(raw), " quantities but ", 1 + n.rep.par + n.nk,
         " were expected; check the ADREPORT lines in ", layout$spec$cpp, ".cpp")
  est <- raw[, "Estimate"]; se <- raw[, "Std. Error"]

  # --- parameters: N = Nmin + exp(Nt); everything else is on the logit scale
  back <- function(x, link) if (link == "N") dat$Nmin + exp(x) else ilogit(x)
  rows <- list(list(labels = "N", link = "N", idx = 1))
  pos  <- 1
  for (b in layout$blocks) {
    rows[[length(rows) + 1]] <- list(labels = b$labels, link = "logit", idx = pos + seq_len(b$n.rep))
    pos <- pos + b$n.rep
  }
  estimates <- do.call(rbind, lapply(rows, function(r)
    data.frame(parameter = r$labels,
               estimate  = back(est[r$idx], r$link),
               lower     = back(est[r$idx] - z * se[r$idx], r$link),
               upper     = back(est[r$idx] + z * se[r$idx], r$link))))
  rownames(estimates) <- NULL

  # --- abundances (count scale, no back-transformation)
  ci <- function(e, s) {
    if (abundance.ci == "wald") cbind(e - z * s, e + z * s)
    else                        cbind(e * exp(-z * s / e), e * exp(z * s / e))
  }
  av  <- pos + seq_len(layout$n.avail)
  al  <- pos + layout$n.avail + seq_len(layout$n.alive)
  tab <- data.frame(period = seq_len(K), Navail = NA_real_, Navail.lo = NA_real_, Navail.up = NA_real_,
                    Nalive = NA_real_, Nalive.lo = NA_real_, Nalive.up = NA_real_)
  tab[, c("Navail", "Navail.lo", "Navail.up")] <- cbind(est[av], ci(est[av], se[av]))
  # alive in period 1 equals available in period 1; later periods come from the C++ code
  alive <- rbind(unlist(tab[1, c("Navail", "Navail.lo", "Navail.up")]),
                 cbind(est[al], ci(est[al], se[al])))
  tab[seq_len(nrow(alive)), c("Nalive", "Nalive.lo", "Nalive.up")] <- alive

  npar <- length(obj$par)
  structure(list(model = model, label = layout$spec$label, penalty = penalty, sigma = sigma,
                 nllike = unname(opt$objective), npar = npar,
                 AIC = 2 * opt$objective + 2 * npar,
                 BIC = 2 * opt$objective + npar * log(dat$Nmin),
                 convergence = opt$convergence, pdHess = isTRUE(sdr$pdHess),
                 estimates = estimates, abundance = tab, raw = raw),
            class = "emig_fit")
}

print.emig_fit <- function(x, digits = 3, ...) {
  cat("\n", x$label, " [", x$model, "]\n", sep = "")
  cat("objective:", round(x$nllike, 3), "| parameters:", x$npar,
      "| AIC:", round(x$AIC, 2), "| BIC:", round(x$BIC, 2),
      "| converged:", x$convergence == 0, "| pd Hessian:", x$pdHess, "\n\n")
  print(format(x$estimates, digits = digits), row.names = FALSE)
  cat("\nAbundance (Navail = available for detection, Nalive = alive):\n")
  print(format(x$abundance, digits = digits), row.names = FALSE)
  invisible(x)
}

#-------------------------------------------------------------------------------
# 9. FIT SEVERAL MODELS AND COMPARE THEM
#-------------------------------------------------------------------------------
# A model that fails is reported and skipped; the others are still fitted.
# save.dir: if given, each model is saved as <save.dir>/<model>.RData (e.g. atbt.RData,
#           containing one object named after the model) as soon as it has been fitted.
fit.models <- function(models = names(MODELS), dat, save.dir = NULL, ...) {
  fits <- list()
  if (!is.null(save.dir)) dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
  for (mod in models) {
    message("Fitting ", mod, " ...")
    res <- tryCatch(fit.model(mod, dat, ...), error = function(e) {
      message("  model '", mod, "' failed: ", conditionMessage(e)); NULL })
    if (!is.null(res)) {
      fits[[mod]] <- res
      if (!is.null(save.dir)) {
        assign(mod, res)
        file <- file.path(save.dir, paste0(mod, ".RData"))
        save(list = mod, file = file)
        message("  saved ", file)
      }
    }
  }
  fits
}

compare.models <- function(fits) {
  if (length(fits) == 0) stop("No fitted models to compare")
  out <- do.call(rbind, lapply(fits, function(f) {
    N <- f$estimates[f$estimates$parameter == "N", ]
    data.frame(model = f$model, objective = f$nllike, npar = f$npar, AIC = f$AIC, BIC = f$BIC,
               N = N$estimate, N.lower = N$lower, N.upper = N$upper,
               converged = f$convergence == 0, pdHess = f$pdHess)
  }))
  rownames(out) <- NULL
  out
}
