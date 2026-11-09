#===============================================================================
# noTE.R  -  the original no-emigration (NoTE) model, phitpt: phi(t), p(t), entry beta(t)
#
# The functions generate.truehists(), find.primhists(), find.event1.indices(),
# find.linkmat() (renamed noTE.linkmat() to avoid a clash with the emigration code)
# and fitLMM.func() are copied unchanged from util_noTE.R. The NoTE model can only
# be fitted to batch-marking data; it is fitted to exactly the same simulated data
# sets as the emigration model.
#===============================================================================


generate.truehists <- function(nOccas){
  ## Return all eligible true capture histories 
  ## nOccas: the length of nOccas is the number of primary periods; 
  ## the k-th element of nOccas is the number of secondary occasions within primary period k
  ## Each true capture history may consist of:
  ## 0 (not captured), 1 (initially captured and marked), and 2 (recaptured)
  ## -------------------------------------------------------------------
  ## Total number of secondary occasions (the length of each true capture history)
  n <- sum(nOccas)
  truehists <- rep(0, n) ## Null history
  ## True histories with event 1 on different occasions
  ## Only 0 before event 1 and only 0 or 2 after event 1
  for(i in 1:n){
    ## Event 1 appears on the i-th occasion
    expandstr <- paste("expand.grid(", paste(c(rep("0", i-1), 1, rep("c(0, 2)", n-i)), collapse=", "), ")")
    truehists <- rbind(truehists, as.matrix(eval(parse(text=expandstr))))
  } 
  colnames(truehists) <- NULL
  return(truehists)
}
##############################################################

find.primhists <- function(truehists, ppInd){
  ## Return corresponding primary histories given true capture histories
  ## For example, given 10 02 00 20 (two occasions for each period), we get 1101
  ## -------------------------------------------------------------------
  ## Find the corresponding primary history for a single true capture history
  find.one.true.primhist <- function(truehist) {
    apply(ppInd, 1, function(ind) ifelse(any(truehist[ind[1]:ind[2]] != 0), 1, 0))
  }
  ## For all histories
  trueprimhists <- t(apply(truehists, 1, find.one.true.primhist))
  return(trueprimhists)
}

#############################################################

find.event1.indices <- function(hists){
  ## Find the indices of the first and last event 1 in each primary history 
  ## For example, given true history=c(1, 0, 1, 0), it will return c(1, 3)
  ## -------------------------------------------------------------------
  inds <- matrix(0, nrow = nrow(hists), ncol = 2)
  for(i in 1:nrow(hists)){
    allinds <- which(hists[i,] == 1)
    if(length(allinds) >= 2) inds[i,] <- c(min(allinds), max(allinds))
    else if(length(allinds) == 1) inds[i,] <- rep(allinds, 2) ## There is only one event 1
    else inds[i,] <- rep(0, 2) ## The null history
  }
  return(inds)
}
##############################################################

noTE.linkmat <- function(truehists, ppInd){
  ## Find the link matrix that connects the latent data to observed data given true capture histories
  ## ppInd[i,] gives starting and ending occasions of primary period i
  ## -------------------------------------------------------------------
  nhist <- nrow(truehists)
  n <- ncol(truehists)   ## Total number of secondary occasions 
  nper <- nrow(ppInd)    ## Total number of primary periods
  
  ## Part I: matrix A (marked individuals)
  Amat <- matrix(0, nrow = n, ncol = nhist)
  for(i in 1:nhist){
    tmp <- which(truehists[i,] == 1)
    if(length(tmp) > 0) Amat[tmp, i] <- 1
  }

  ## Part II: matrix B (recaptured individuals)
  Bmat <- matrix(0, nrow = nper * n, ncol = nhist)
  for(i in 1:nhist){
    ## Check if there are any recaptures (event 2) first
    recap <- any(truehists[i,] == 2)
    if(recap){
      ## The occasion on which the individual is marked
      markocca <- which(truehists[i,] == 1)
      tmp <- ppInd - markocca 
      ## The corresponding primary period
      markpp <- which(tmp[,1]*tmp[,2] <= 0)
      indstart <- (markpp - 1) * n + markocca + 1
      indend <- markpp * n 
      for(j in indstart:indend){
        if(truehists[i, j-n*(markpp-1)] == 2) Bmat[j, i] <- 1
      }
    }
  }
  ## Some rows of Cmat always comprise of zero entries and thus are not needed in the framework
  zeroRowInd <- which(rowSums(Bmat) == 0)
  Bmat <- Bmat[-zeroRowInd, ]
  
  ## Return the combination of the three matrices
  linkmat <- rbind(Amat, Bmat)
  return(linkmat)
}
#############################################################

fitLMM.func <- function(obscounts, linkmat, truehists, trueprimhists, ppInd, eventoneInd, 
                        startvec = NULL, model = "phitpt", penalty = 1, sigma = 1.25){
  ## Fit the latent multinomial model to extended batach mark data 
  ## startvec is a vector of starting values for optimization
  ## startvec = c(N.start, p.start, phi.start, beta0.start)
  ## ------------------------------------------------------------------- 
  nocca <- ncol(truehists) ## Total number of capture occasions
  nper  <- ncol(trueprimhists) ## Number of primary periods
  ## Check the inputs
  # if(length(obscounts) != nrow(linkmat)       | 
  #    ncol(linkmat)     != nrow(truehists)     | 
  #    nrow(truehists)   != nrow(trueprimhists) |
  #    nrow(truehists)   != nrow(eventoneInd)   |
  #    nrow(ppInd)       != nper) stop("Incorrect inputs detected")
  ## Determine a lower bound for N
  ## N should not be smaller than the total number of individuals marked
  Nmin <- sum(obscounts[1:nocca])  
  
  ## If there are zero counts in the observed vector, some adjustment is needed to ensure that the SPA works
  ## This will reduce the amount of computation needed and improve efficiency
  if(any(obscounts == 0)){
    zero.counts.ind <- which(obscounts == 0)
    zero.truehists.ind <- which(colSums(linkmat[zero.counts.ind,,drop = F]) > 0)
    linkmat <- linkmat[-zero.counts.ind, -zero.truehists.ind]
    truehists <- truehists[-zero.truehists.ind,]
    trueprimhists <- trueprimhists[-zero.truehists.ind,]
    eventoneInd <- eventoneInd[-zero.truehists.ind,]
    obscounts <- obscounts[-zero.counts.ind]
  }
  nhist <- nrow(truehists) ## Number of true capture histories
  ## Consider different models
  if(model == "phitpt"){
    np <- nocca; nphi <- nper -1 
  }else if(model == "phitpc"){
    np <- 1; nphi <- nper -1 
  }else if(model == "phicpt"){
    np <- nocca; nphi <- 1 
  }else if(model == "phicpc"){
    np <- 1; nphi <- 1 
  }else {
    stop("A wrong model name provided.")
  }
  ## Sort out start values for optimisation
  if(is.null(startvec)){
    Nt.start     <- 1
    pt.start     <- rep(1, np)
    phit.start   <- rep(1, nphi)
    beta0t.start <- rep(1, nper-1)
  }else {
    Nt.start     <- log(startvec[1] - Nmin)
    pt.start     <- logit(startvec[2:(np+1)])
    phit.start   <- logit(startvec[(np+2):(nphi+np+1)])
    beta0t.start <- trans.beta.func(startvec[(nphi+np+2):(np+nphi+nper)])
  }
  ## Penalty items
  phi_lambda   <- rep(1/(2 * sigma^2), nphi)
  p_lambda     <- rep(1/(2 * sigma^2), np)
  beta0_lambda <- rep(1/(2 * sigma^2), (nper-1))

  ## Construct AD template function
  dat <- list(obs_counts = obscounts, linkmat = as(linkmat, "dgTMatrix"),
              true_hists = truehists, true_prim_hists = trueprimhists, 
              pp_ind = ppInd, eventone_ind = eventoneInd, Nmin = Nmin,
              nocca = nocca, nper = nper, nhist = nhist,
              phi_lambda = phi_lambda, p_lambda = p_lambda, beta0_lambda = beta0_lambda,
              penalty = penalty)
  pars <- list(Nt = Nt.start, pt = pt.start, phit = phit.start, 
               beta0t = beta0t.start, svec = rep(0, length(obscounts)))
  obj <- MakeADFun(dat, pars, random = "svec", DLL = model, silent = TRUE,
                   saddlepoint = TRUE)
  ## Fit the model using the optimizer nlminb
  catline("Start fitting", paste("model", model), "...")
  time <- system.time(optres <- nlminb(obj$par, obj$fn, obj$gr,
                                       control = list(iter.max=2000, eval.max=5000)))
  catline("Done!")
  catline("Time used for optimization:", paste(round(time[1], 2), "s", sep=""))
  ## Convergence code: 0 means successful convergence
  convg <- optres$convergence
  ## If the convergence code is non-zero, do one more optimization
  if(convg != 0)  optres <- nlminb(optres$par, obj$fn, obj$gr,
                                   control = list(iter.max=2000, eval.max=5000))
  convg <- optres$convergence
  ## If it still does not converge, throw out a warning message
  if(convg !=0 ) warning("Convergence code of nlminb is non-zero")
  names(convg) <- "convg"
  
  ## Negative loglikelihood
  nllike <- optres$objective
  names(nllike) <- "nllike"
  
  ## Find variances of the estimates.
  ## Note that calling sdreport(obj) now returns the sd report at the MLEs.
  est.details <- sdreport(obj)
  est.mat <- summary.sdreport(est.details, select = c("report"))
  
  # extract primary period abundance estimates 
  nk_ind <- length(est.mat[,"Estimate"])-(nper)+1
  nk_est <- est.mat[,"Estimate"][nk_ind:(length(est.mat[,"Estimate"]))]
  nk_se <- est.mat[,"Std. Error"][nk_ind:(length(est.mat[,"Estimate"]))]
  nk_Lo <- nk_est - 1.959964 * nk_se
  nk_Up <- nk_est + 1.959964 * nk_se 
  
  transMLEs <- est.mat[,"Estimate"][-c(nk_ind:(length(est.mat[,"Estimate"])))]
  ## Std. Errors
  se.ests <- est.mat[,"Std. Error"][-c(nk_ind:(length(est.mat[,"Estimate"])))]
  
  if(any(is.nan(se.ests))) warning("TMB returns NaN standard errors")
  ## 95% CIs
  transCI.Lo <- transMLEs - 1.959964 * se.ests
  transCI.Up <- transMLEs + 1.959964 * se.ests 
  ## Back transform the parameters and CIs
  N.est <- Nmin + exp(transMLEs[1])
  N.Lo  <- Nmin + exp(transCI.Lo[1])
  N.Up  <- Nmin + exp(transCI.Up[1])
  ## Other parameters
  param.est <- ilogit(transMLEs[-1])
  param.Lo  <- ilogit(transCI.Lo[-1])
  param.Up  <- ilogit(transCI.Up[-1])
  ## Put MLEs and CIs together
  MLEs <- c(N.est, param.est)
  CI.Lo <- c(N.Lo, param.Lo)
  CI.Up <- c(N.Up, param.Up)
  ## Name the output
  names(MLEs) <- c("N", 
                   paste0("p", 1:np), 
                   paste0("phi", 1:nphi),
                   paste0("beta", 1:nper))
  names(CI.Lo) <- paste0(names(MLEs), ".Lo")
  names(CI.Up) <- paste0(names(MLEs), ".Up")
  ## Output
  allres <- c(MLEs, nk_est, CI.Lo, nk_Lo, CI.Up, nk_Up, nllike, convg)
  ## Note that exp(x)/(1+exp(x)) is NaN if x is large. 
  ## This should be 1 instead.
  allres[which(is.nan(allres))] <- 1 
  return(allres)
}


#===============================================================================
# Wrappers used by the simulation scripts
#===============================================================================

# NoTE latent structure for a design (slow for many occasions, so it is cached)
noTE.structure <- function(secondary, cache.dir = NULL) {
  file <- if (!is.null(cache.dir))
    file.path(cache.dir, paste0("latent_", paste(secondary, collapse = "-"), "_noTE.rds"))
  if (!is.null(file) && file.exists(file)) {
    message("Reading NoTE latent structure from ", file)
    return(readRDS(file))
  }
  message("Generating NoTE latent histories and link matrix (slow, done only once)...")
  ppInd         <- find.primperiod.indices(secondary)
  truehists     <- generate.truehists(secondary)
  trueprimhists <- find.primhists(truehists, ppInd)
  str <- list(secondary = secondary, ppInd = ppInd, truehists = truehists,
              trueprimhists = trueprimhists,
              eventoneInd = find.event1.indices(trueprimhists),
              linkmat = noTE.linkmat(truehists, ppInd))
  if (!is.null(file)) {
    dir.create(cache.dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(str, file)
  }
  str
}

# Fit the NoTE model (phitpt) to ONE simulated data set and return a one-row data.frame
# in the same layout as the emigration results: iteration, objective, convergence,
# estimates, lower limits (.Lo), upper limits (.Up).
fit.noTE.once <- function(obscounts, iteration, str, S, penalty = 0) {

  nocca <- sum(S$n.secondary)
  nper  <- length(S$n.secondary)
  if (length(obscounts) != nrow(str$linkmat))
    stop("The observed data have ", length(obscounts), " counts but the NoTE link matrix has ",
         nrow(str$linkmat), " rows")

  start <- c(S$N, rep(S$p, nocca), rep(S$phi, nper - 1), rep(1 / nper, nper - 1))
  res <- fitLMM.func(obscounts, str$linkmat, str$truehists, str$trueprimhists, str$ppInd,
                     str$eventoneInd, start, model = "phitpt", penalty = penalty, sigma = S$sigma)

  nMLE  <- 1 + nocca + (nper - 1) + nper
  names.par <- c("N", paste0("p", seq_len(nocca)), paste0("phi", seq_len(nper - 1)),
                 paste0("beta", seq_len(nper)))
  names.N   <- paste0("Nalive", seq_len(nper))
  nm  <- c(names.par, names.N)
  # fitLMM.func returns c(MLEs, N_alive, lower(MLEs), lower(N_alive), upper(MLEs), upper(N_alive), nllike, convg)
  blk <- nMLE + nper
  est <- res[seq_len(blk)]
  lo  <- res[blk + seq_len(blk)]
  up  <- res[2 * blk + seq_len(blk)]
  out <- c(iteration = iteration, objective = unname(res[3 * blk + 1]),
           convergence = unname(res[3 * blk + 2]),
           setNames(est, nm), setNames(lo, paste0(nm, ".Lo")), setNames(up, paste0(nm, ".Up")))
  as.data.frame(as.list(out), check.names = FALSE)
}
