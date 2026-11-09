#===============================================================================
# emigration_models.R  -  the five emigration models: registry, parameter layout,
#                         compile/load, and a function that fits ONE model to ONE
#                         data set (re-used from the real-data analysis).
#
# Models (the name is what you pass to the functions):
#   "cr"   completely random emigration, time-dependent alpha
#   "atbt" Markovian emigration, time-dependent alpha, time-dependent beta
#   "atbc" Markovian emigration, time-dependent alpha, constant beta
#   "acbt" Markovian emigration, constant alpha,       time-dependent beta
#   "acbc" Markovian emigration, constant alpha,       constant beta
#===============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(TMB)
})

#-------------------------------------------------------------------------------
# MODEL REGISTRY
#    "confounded" models are those where the last survival/emigration (and return)
#    parameters cannot be separated, so only the products phi*(1-alpha) (and phi*beta)
#    are estimated.
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
# DATA PREPARATION
#-------------------------------------------------------------------------------
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
# MODEL LAYOUT: which parameters each model has, how many, and their starting values.
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
# COMPILE / LOAD THE C++ FILES
#-------------------------------------------------------------------------------
# cpp.names: names of the C++ files without extension, e.g. "markovian_atbt" or "phitpt".
# Compiles only if the compiled file is missing or older than the source, then loads it.
load.cpp <- function(cpp.names, dir = "cpp", flags = "-O1 -g") {
  old <- setwd(dir); on.exit(setwd(old))
  for (cpp in cpp.names) {
    src <- paste0(cpp, ".cpp")
    if (!file.exists(src)) stop("Cannot find ", file.path(dir, src))
    if (!file.exists(dynlib(cpp)) || file.mtime(dynlib(cpp)) < file.mtime(src))
      compile(src, flags, DLLFLAGS = "")
    if (!cpp %in% names(getLoadedDLLs())) dyn.load(dynlib(cpp))
  }
  invisible(cpp.names)
}

#-------------------------------------------------------------------------------
# FIT ONE MODEL
#-------------------------------------------------------------------------------
# dat       list with nprim, nocca, ppInd, true_hists, linkmat, obscounts, Nmin
# penalty   1 = penalised likelihood (ridge penalty on the logit parameters), 0 = none
# sigma     penalty strength: lambda = 1 / (2 sigma^2)
# start     starting values (probabilities) for p, phi, alpha, beta
# Nt.start  starting value for log(N - Nmin)
# abundance.ci  "wald" (symmetric, count scale) or "log" (log-normal) intervals for Navail / Nalive
fit.model <- function(model, dat, penalty = 1, sigma = 2.5, start = default.start(),
                      Nt.start = 1, abundance.ci = c("wald", "log"), silent = TRUE) {

  abundance.ci <- match.arg(abundance.ci)
  z       <- Z.CRIT
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
