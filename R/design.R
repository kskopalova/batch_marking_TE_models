#===============================================================================
# design.R  -  latent histories and link matrices (depend on the design only)
#
#   Latent states: 0 not yet entered, 1 present/not captured, 2 first capture,
#                  3 temporarily away, 4 dead  (5 = recapture, used in the link matrix)
#
#   data.type = "batch"       batch-marking: observed data are c(m, n)
#                             (first captures per occasion, recaptures per marking period)
#   data.type = "individual"  individual identification: observed data are the counts
#                             of every observed capture history
#===============================================================================

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

# Batch-marking link matrix: rows = observed counts c(m, n), columns = latent histories.
# (Same matrix as linkmat.batchmark() in the original util_emig.R.)
linkmat.batchmark <- function(th, sc) {
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


# Individual-identification link matrix: one row per possible observed capture history
# (rows ordered as expand.grid(rep(list(0:1), nocca)) without the all-zero history),
# one column per latent history; an entry is 1 if the latent history is observed as that
# capture history. Same matrix as linkmat.ind() in the original util_emig.R, built directly
# from the binary code of the capture history instead of comparing all pairs.
linkmat.individual <- function(th) {
  nocca <- ncol(th)
  code  <- as.vector((th == 2) %*% (2^(seq_len(nocca) - 1)))   # row index of the observed history
  seen  <- which(code > 0)                                     # never-captured histories are not observed
  sparseMatrix(i = code[seen], j = seen, x = rep(1, length(seen)),
               dims = c(2^nocca - 1, nrow(th)))
}

# Latent histories + link matrix for a design and data type. They are slow to build and only
# depend on the design, so they are stored in cache.dir and re-used on later runs.
latent.structure <- function(secondary, data.type = c("batch", "individual"), cache.dir = NULL) {
  data.type <- match.arg(data.type)
  file <- if (!is.null(cache.dir))
    file.path(cache.dir, paste0("latent_", paste(secondary, collapse = "-"), "_", data.type, ".rds"))
  if (!is.null(file) && file.exists(file)) {
    message("Reading latent histories and link matrix from ", file)
    return(readRDS(file))
  }
  message("Generating latent histories and link matrix (slow, done only once)...")
  th  <- generate.latent.histories(secondary)
  lat <- list(secondary = secondary, data.type = data.type, true_hists = th,
              linkmat = if (data.type == "batch") linkmat.batchmark(th, secondary) else linkmat.individual(th))
  if (!is.null(file)) {
    dir.create(cache.dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(lat, file)
  }
  lat
}
