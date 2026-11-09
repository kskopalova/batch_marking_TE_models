#===============================================================================
# true_model.R
#
# The 'true' data-generating models: latent-history probabilities for each
# emigration model and the functions that turn the simulation settings into
# parameter vectors (xxx.assign) and latent probabilities (xxx.probs).
#
# Everything in this file is copied unchanged from the original util_emig.R
# (transition matrices, latent.probs.calc.*, *.assign, *.probs); only the
# unused runsims.* / analyse.* / run.full.sim.study.* functions were left out.
#===============================================================================

within.primary.transitions <- function(p) {
  
  x <- matrix(nrow=5, ncol=5)
  rownames(x) <- paste0("state", 0:4)
  colnames(x) <- paste0("state", 0:4)
  
  x[1,] <- c(1, 0, 0, 0, 0)
  x[2,] <- c(0, 1-p, p, 0, 0)
  x[3,] <- c(0, 1-p, p, 0, 0)
  x[4,] <- c(0, 0, 0, 1, 0)
  x[5,] <- c(0, 0, 0, 0, 1)
  return(as.matrix(x))
}

between.primary.transitions <- function(mig.prob, 
                                        return.prob, 
                                        entry.prob, 
                                        surv.prob, 
                                        capture.prob) {
  
  x <- matrix(nrow=5, ncol=5)
  rownames(x) <- paste0("state", 0:4)
  colnames(x) <- paste0("state", 0:4)
  
  x[1,] <- c(1, entry.prob*(1-capture.prob), entry.prob*capture.prob, 0, 0)
  x[2,] <- c(0, surv.prob*(1-mig.prob)*(1-capture.prob), surv.prob*(1-mig.prob)*capture.prob, surv.prob*mig.prob, 1-surv.prob)
  x[3,] <- c(0, surv.prob*(1-mig.prob)*(1-capture.prob), surv.prob*(1-mig.prob)*capture.prob, surv.prob*mig.prob, 1-surv.prob)
  x[4,] <- c(0, surv.prob*return.prob*(1-capture.prob), surv.prob*return.prob*capture.prob, surv.prob*(1-return.prob), 1-surv.prob)
  x[5,] <- c(rep(0, 4), 1)
  return(as.matrix(x))
}

last.two.periods <- function(mig.prob, 
                             return.prob, 
                             entry.prob, 
                             surv.prob, 
                             capture.prob,
                             mig.prod,
                             return.prod) {
  
  x <- matrix(nrow=5, ncol=5)
  rownames(x) <- paste0("state", 0:4)
  colnames(x) <- paste0("state", 0:4)
  
  x[1,] <- c(1, entry.prob*(1-capture.prob), entry.prob*capture.prob, 0, 0)
  x[2,] <- c(0, mig.prod*(1-capture.prob), mig.prod*capture.prob, (1-mig.prod)/2, (1-mig.prod)/2)
  x[3,] <- c(0, mig.prod*(1-capture.prob), mig.prod*capture.prob, (1-mig.prod)/2, (1-mig.prod)/2)
  x[4,] <- c(0, return.prod*(1-capture.prob), return.prod*capture.prob, (1-return.prod)/2, (1-return.prod)/2)
  x[5,] <- c(rep(0, 4), 1)
  return(as.matrix(x))
}

# completely random emigration model:
# uses the same within-transition matrix
between.primary.transitions.cr <- function(mig.prob,  
                                           entry.prob, 
                                           surv.prob, 
                                           capture.prob) {
  
  x <- matrix(nrow=5, ncol=5)
  rownames(x) <- paste0("state", 0:4)
  colnames(x) <- paste0("state", 0:4)
  
  x[1,] <- c(1, entry.prob*(1-capture.prob), entry.prob*capture.prob, 0, 0)
  x[2,] <- c(0, surv.prob*(1-mig.prob)*(1-capture.prob), surv.prob*(1-mig.prob)*capture.prob, surv.prob*mig.prob, 1-surv.prob)
  x[3,] <- c(0, surv.prob*(1-mig.prob)*(1-capture.prob), surv.prob*(1-mig.prob)*capture.prob, surv.prob*mig.prob, 1-surv.prob)
  x[4,] <- c(0, surv.prob*(1-mig.prob)*(1-capture.prob), surv.prob*(1-mig.prob)*capture.prob, surv.prob*mig.prob, 1-surv.prob)
  x[5,] <- c(rep(0, 4), 1)
  return(as.matrix(x))
}

last.two.periods.cr <- function(mig.prob, 
                                entry.prob, 
                                surv.prob, 
                                capture.prob,
                                mig.prod) {
  
  x <- matrix(nrow=5, ncol=5)
  rownames(x) <- paste0("state", 0:4)
  colnames(x) <- paste0("state", 0:4)
  
  x[1,] <- c(1, entry.prob*(1-capture.prob), entry.prob*capture.prob, 0, 0)
  x[2,] <- c(0, mig.prod*(1-capture.prob), mig.prod*capture.prob, (1-mig.prod)/2, (1-mig.prod)/2)
  x[3,] <- c(0, mig.prod*(1-capture.prob), mig.prod*capture.prob, (1-mig.prod)/2, (1-mig.prod)/2)
  x[4,] <- c(0, mig.prod*(1-capture.prob), mig.prod*capture.prob, (1-mig.prod)/2, (1-mig.prod)/2)
  x[5,] <- c(rep(0, 4), 1)
  return(as.matrix(x))
}


latent.probs.calc.atbt <- function(nocca, true_hists, 
                              p, phi, gamma, alpha, beta, prodalpha, prodbeta,
                              between.primary.ind) {
  
  lst <- vector(mode = "list", length = nocca-1)
  names(lst) <- paste0("occasion", 1:(nocca-1), "to", 2:nocca)
  
  # assign transition matrices to the list
  for (i in 1:length(lst)) {
    
    if (i==between.primary.ind[length(between.primary.ind)-1]) {
      lst[[i]] <- last.two.periods(mig.prob=alpha[i], 
                                   return.prob=beta[i], 
                                   entry.prob=gamma[i+1],
                                   surv.prob=phi[i], 
                                   capture.prob=p[i+1],
                                   mig.prod=prodalpha,
                                   return.prod=prodbeta)
    } else if (i %in% between.primary.ind) {
      lst[[i]] <- between.primary.transitions(mig.prob=alpha[i], 
                                              return.prob=beta[i], 
                                              entry.prob=gamma[i+1],
                                              surv.prob=phi[i], 
                                              capture.prob=p[i+1])
    } else {
      lst[[i]] <- within.primary.transitions(p=p[i+1])
    }
  }
  
  arrival.ind <- apply(true_hists, 1, function(x) which(x==1 | x==2)[1])
  arrival.state <- rep(NA, nrow(true_hists))
  for (i in 1:length(arrival.state)) {
    arrival.state[i] <- ifelse(is.na(arrival.ind[i]), NA, as.numeric(true_hists[i,][arrival.ind[i]]))
  }
  
  true_hists <- true_hists[!is.na(arrival.ind),]
  probs <- rep(NA, nrow(true_hists))
  arrival.ind <- arrival.ind[!is.na(arrival.ind)]
  arrival.state <- arrival.state[!is.na(arrival.state)]
  
  # this counts from occasion 1
  probs <- rep(1, nrow(true_hists))
  for (i in 1:nrow(true_hists)) {
    for (j in 1:(nocca-1)) {
      state.t <- as.numeric(true_hists[i,][j])
      state.tplusone <- as.numeric(true_hists[i,][j+1])
      
      tmp <- lst[[j]][rownames(lst[[j]])==paste0("state", state.t),
                      colnames(lst[[j]])==paste0("state", state.tplusone)]
      probs[i] <- probs[i]*tmp
    }
  }
  for (i in 1:length(probs)) {
    
    if (true_hists[i,1] == 1 | true_hists[i,1] == 2) {
      
      prob.at.arrival <- ifelse(arrival.state[i]==1, 
                                gamma[1]*(1-p[1]), 
                                gamma[1]*p[1])
      probs[i] <- probs[i]*prob.at.arrival
    } else if (true_hists[i,1] == 4 | true_hists[i,1] == 3) {
      probs[i] <- 0
    }
    
  }
  return(probs)
}

latent.probs.calc.atbc <- function(true_hists, nocca, p, phi, gamma, alpha, beta, between.primary.ind) {
  
  lst <- vector(mode = "list", length = nocca-1)
  names(lst) <- paste0("occasion", 1:(nocca-1), "to", 2:nocca)
  
  # assign transition matrices to the list
  for (i in 1:length(lst)) {
    
    if (i %in% between.primary.ind) {
      lst[[i]] <- between.primary.transitions(mig.prob=alpha[i], 
                                              return.prob=beta[i], 
                                              entry.prob=gamma[i+1],
                                              surv.prob=phi[i], 
                                              capture.prob=p[i+1])
    } else {
      lst[[i]] <- within.primary.transitions(p=p[i+1])
    }
  }
  
  arrival.ind <- apply(true_hists, 1, function(x) which(x==1 | x==2)[1])
  arrival.state <- rep(NA, nrow(true_hists))
  for (i in 1:length(arrival.state)) {
    arrival.state[i] <- ifelse(is.na(arrival.ind[i]), NA, as.numeric(true_hists[i,][arrival.ind[i]]))
  }
  
  true_hists <- true_hists[!is.na(arrival.ind),]
  probs <- rep(NA, nrow(true_hists))
  arrival.ind <- arrival.ind[!is.na(arrival.ind)]
  arrival.state <- arrival.state[!is.na(arrival.state)]
  
  # this counts from occasion 1
  probs <- rep(1, nrow(true_hists))
  for (i in 1:nrow(true_hists)) {
    for (j in 1:(nocca-1)) {
      state.t <- as.numeric(true_hists[i,][j])
      state.tplusone <- as.numeric(true_hists[i,][j+1])
      
      tmp <- lst[[j]][rownames(lst[[j]])==paste0("state", state.t),
                      colnames(lst[[j]])==paste0("state", state.tplusone)]
      probs[i] <- probs[i]*tmp
    }
  }
  for (i in 1:length(probs)) {
    
    if (true_hists[i,1] == 1 | true_hists[i,1] == 2) {
      
      prob.at.arrival <- ifelse(arrival.state[i]==1, 
                                gamma[1]*(1-p[1]), 
                                gamma[1]*p[1])
      probs[i] <- probs[i]*prob.at.arrival
    } else if (true_hists[i,1] == 4 | true_hists[i,1] == 3) {
      probs[i] <- 0
    }
    
  }
  return(probs)
}

latent.probs.calc.acbt <- function(true_hists, nocca, p, phi, gamma, alpha, beta, between.primary.ind) {
  
  lst <- vector(mode = "list", length = nocca-1)
  names(lst) <- paste0("occasion", 1:(nocca-1), "to", 2:nocca)
  
  # assign transition matrices to the list
  for (i in 1:length(lst)) {
    
    if (i %in% between.primary.ind) {
      lst[[i]] <- between.primary.transitions(mig.prob=alpha[i], 
                                              return.prob=beta[i], 
                                              entry.prob=gamma[i+1],
                                              surv.prob=phi[i], 
                                              capture.prob=p[i+1])
    } else {
      lst[[i]] <- within.primary.transitions(p=p[i+1])
    }
  }
  
  arrival.ind <- apply(true_hists, 1, function(x) which(x==1 | x==2)[1])
  arrival.state <- rep(NA, nrow(true_hists))
  for (i in 1:length(arrival.state)) {
    arrival.state[i] <- ifelse(is.na(arrival.ind[i]), NA, as.numeric(true_hists[i,][arrival.ind[i]]))
  }
  
  true_hists <- true_hists[!is.na(arrival.ind),]
  probs <- rep(NA, nrow(true_hists))
  arrival.ind <- arrival.ind[!is.na(arrival.ind)]
  arrival.state <- arrival.state[!is.na(arrival.state)]
  
  # this counts from occasion 1
  probs <- rep(1, nrow(true_hists))
  for (i in 1:nrow(true_hists)) {
    for (j in 1:(nocca-1)) {
      state.t <- as.numeric(true_hists[i,][j])
      state.tplusone <- as.numeric(true_hists[i,][j+1])
      
      tmp <- lst[[j]][rownames(lst[[j]])==paste0("state", state.t),
                      colnames(lst[[j]])==paste0("state", state.tplusone)]
      probs[i] <- probs[i]*tmp
    }
  }
  for (i in 1:length(probs)) {
    
    if (true_hists[i,1] == 1 | true_hists[i,1] == 2) {
      
      prob.at.arrival <- ifelse(arrival.state[i]==1, 
                                gamma[1]*(1-p[1]), 
                                gamma[1]*p[1])
      probs[i] <- probs[i]*prob.at.arrival
    } else if (true_hists[i,1] == 4 | true_hists[i,1] == 3) {
      probs[i] <- 0
    }
    
  }
  return(probs)
}

latent.probs.calc.acbc <- function(true_hists, nocca, p, phi, gamma, alpha, beta, between.primary.ind) {
  
  lst <- vector(mode = "list", length = nocca-1)
  names(lst) <- paste0("occasion", 1:(nocca-1), "to", 2:nocca)
  
  # assign transition matrices to the list
  for (i in 1:length(lst)) {
    
    if (i %in% between.primary.ind) {
      lst[[i]] <- between.primary.transitions(mig.prob=alpha[i], 
                                              return.prob=beta[i], 
                                              entry.prob=gamma[i+1],
                                              surv.prob=phi[i], 
                                              capture.prob=p[i+1])
    } else {
      lst[[i]] <- within.primary.transitions(p=p[i+1])
    }
  }
  
  arrival.ind <- apply(true_hists, 1, function(x) which(x==1 | x==2)[1])
  arrival.state <- rep(NA, nrow(true_hists))
  for (i in 1:length(arrival.state)) {
    arrival.state[i] <- ifelse(is.na(arrival.ind[i]), NA, as.numeric(true_hists[i,][arrival.ind[i]]))
  }
  
  true_hists <- true_hists[!is.na(arrival.ind),]
  probs <- rep(NA, nrow(true_hists))
  arrival.ind <- arrival.ind[!is.na(arrival.ind)]
  arrival.state <- arrival.state[!is.na(arrival.state)]
  
  # this counts from occasion 1
  probs <- rep(1, nrow(true_hists))
  for (i in 1:nrow(true_hists)) {
    for (j in 1:(nocca-1)) {
      state.t <- as.numeric(true_hists[i,][j])
      state.tplusone <- as.numeric(true_hists[i,][j+1])
      
      tmp <- lst[[j]][rownames(lst[[j]])==paste0("state", state.t),
                      colnames(lst[[j]])==paste0("state", state.tplusone)]
      probs[i] <- probs[i]*tmp
    }
  }
  for (i in 1:length(probs)) {
    
    if (true_hists[i,1] == 1 | true_hists[i,1] == 2) {
      
      prob.at.arrival <- ifelse(arrival.state[i]==1, 
                                gamma[1]*(1-p[1]), 
                                gamma[1]*p[1])
      probs[i] <- probs[i]*prob.at.arrival
    } else if (true_hists[i,1] == 4 | true_hists[i,1] == 3) {
      probs[i] <- 0
    }
    
  }
  return(probs)
}

latent.probs.calc.cr <- function(nocca, true_hists, 
                              p, phi, gamma, alpha, prodalpha,
                              between.primary.ind) {
  
  lst <- vector(mode = "list", length = nocca-1)
  names(lst) <- paste0("occasion", 1:(nocca-1), "to", 2:nocca)
  
  # assign transition matrices to the list
  for (i in 1:length(lst)) {
    
    if (i==between.primary.ind[length(between.primary.ind)-1]) {
      lst[[i]] <- last.two.periods.cr(mig.prob=alpha[i], 
                                   entry.prob=gamma[i+1],
                                   surv.prob=phi[i], 
                                   capture.prob=p[i+1],
                                   mig.prod=prodalpha)
    } else if (i %in% between.primary.ind) {
      lst[[i]] <- between.primary.transitions.cr(mig.prob=alpha[i], 
                                              entry.prob=gamma[i+1],
                                              surv.prob=phi[i], 
                                              capture.prob=p[i+1])
    } else {
      lst[[i]] <- within.primary.transitions(p=p[i+1])
    }
  }
  
  arrival.ind <- apply(true_hists, 1, function(x) which(x==1 | x==2)[1])
  arrival.state <- rep(NA, nrow(true_hists))
  for (i in 1:length(arrival.state)) {
    arrival.state[i] <- ifelse(is.na(arrival.ind[i]), NA, as.numeric(true_hists[i,][arrival.ind[i]]))
  }
  
  true_hists <- true_hists[!is.na(arrival.ind),]
  probs <- rep(NA, nrow(true_hists))
  arrival.ind <- arrival.ind[!is.na(arrival.ind)]
  arrival.state <- arrival.state[!is.na(arrival.state)]
  
  # this counts from occasion 1
  probs <- rep(1, nrow(true_hists))
  for (i in 1:nrow(true_hists)) {
    for (j in 1:(nocca-1)) {
      state.t <- as.numeric(true_hists[i,][j])
      state.tplusone <- as.numeric(true_hists[i,][j+1])
      
      tmp <- lst[[j]][rownames(lst[[j]])==paste0("state", state.t),
                      colnames(lst[[j]])==paste0("state", state.tplusone)]
      probs[i] <- probs[i]*tmp
    }
  }
  for (i in 1:length(probs)) {
    
    if (true_hists[i,1] == 1 | true_hists[i,1] == 2) {
      
      prob.at.arrival <- ifelse(arrival.state[i]==1, 
                                gamma[1]*(1-p[1]), 
                                gamma[1]*p[1])
      probs[i] <- probs[i]*prob.at.arrival
    } else if (true_hists[i,1] == 4 | true_hists[i,1] == 3) {
      probs[i] <- 0
    }
    
  }
  return(probs)
}

#----- model-specific parameter assignment and latent probabilities --------

atbt.assign <- function(seed, 
                        prim, 
                        occa, 
                        between.inds, 
                        lst_param=vector(mode="list", length=7), 
                        lst_start=vector(mode="list", length=7),
                        migrationprob,
                        returnprob,
                        captureprob,
                        survivalprob) {
  
  set.seed(seed)
  
  aprim <- rep(migrationprob, prim-1)
  last.alpha <- aprim[length(aprim)]
  aprim <- aprim[-length(aprim)]
  a <- rep(1, occa-1)
  a[cumsum(between.inds)[1:length(aprim)]] <- aprim
  
  bprim <- rep(returnprob, prim-2)
  last.beta <- bprim[length(bprim)]
  bprim <- bprim[-length(bprim)]
  b <- rep(1, occa-1)
  b[cumsum(between.inds)[-c(1, prim-1, prim)]] <- bprim
  
  p <- rep(captureprob, occa)
  
  phiprim <- rep(survivalprob, prim-1)
  last.phi <- phiprim[length(phiprim)]
  phiprim <- phiprim[-length(phiprim)]
  phi <- rep(1, occa-1)
  phi[cumsum(between.inds)[-c(prim-1, prim)]] <- phiprim
  
  gprim <- rep(1/(prim), prim)
  g <- rep(1, occa)
  g[c(1,cumsum(between.inds)+1)[-prim-1]] <- gprim
  
  # prod.alpha <- (1-last.phi*(1-last.alpha))/2
  # prod.beta <- (1-last.phi*last.beta)/2
  
  prod.alpha <- last.phi*(1-last.alpha)
  prod.beta <- last.phi*last.beta
  
  lst_param[[1]] <- a
  lst_param[[2]] <- b
  lst_param[[3]] <- p
  lst_param[[4]] <- phi
  lst_param[[5]] <- g
  lst_param[[6]] <- prod.alpha
  lst_param[[7]] <- prod.beta
  
  lst_start[[1]] <- logit(aprim)
  lst_start[[2]] <- logit(bprim)
  lst_start[[3]] <- logit(p)
  lst_start[[4]] <- logit(phiprim)
  g0 <- gprim[-prim] 
  lst_start[[5]] <- trans.beta.func(g0)
  lst_start[[6]] <- logit(prod.alpha)
  lst_start[[7]] <- logit(prod.beta)
  
  return(list(atbt_param=lst_param, atbt_starting=lst_start))
  
}


atbt.probs <- function(th, param_list, nOccas, sc) {
  
  probs <- rep(NA, nrow(th))
  for (i in 1:nrow(th)) {
    probs[i] <- latent.probs.calc.atbt(nocca=nOccas,
                                       true_hists=th[i,], 
                                       p=param_list[[3]],
                                       phi=param_list[[4]],
                                       gamma=param_list[[5]],
                                       alpha=param_list[[1]],
                                       beta=param_list[[2]],
                                       prodalpha=param_list[[6]],
                                       prodbeta=param_list[[7]],
                                       between.primary.ind=cumsum(sc))
  }
  return(probs)
}

atbc.assign <- function(seed, 
                        prim, 
                        occa, 
                        between.inds, 
                        lst_param=vector(mode="list", length=5), 
                        lst_start=vector(mode="list", length=5),
                        migrationprob,
                        returnprob,
                        captureprob,
                        survivalprob) {
  
  set.seed(seed)
  
  aprim <- rep(migrationprob, prim-1)
  a <- rep(1, occa-1)
  a[cumsum(between.inds)[-prim]] <- aprim
  
  bprim <- returnprob
  b <- rep(1, occa-1)
  b[cumsum(between.inds)[-c(1,prim)]] <- bprim
  
  p <- rep(captureprob, occa)
  
  phiprim <- rep(survivalprob, prim-1)
  phi <- rep(1, occa-1)
  phi[cumsum(between.inds)[-prim]] <- phiprim
  
  gprim <- rep(1/(prim), prim)
  g <- rep(1, occa)
  g[c(1,cumsum(between.inds)+1)[-prim-1]] <- gprim
  
  lst_param[[1]] <- a
  lst_param[[2]] <- b
  lst_param[[3]] <- p
  lst_param[[4]] <- phi
  lst_param[[5]] <- g
  
  lst_start[[1]] <- logit(aprim)
  lst_start[[2]] <- logit(bprim)
  lst_start[[3]] <- logit(p)
  lst_start[[4]] <- logit(phiprim)
  g0 <- gprim[-prim] 
  lst_start[[5]] <- trans.beta.func(g0)
  
  return(list(atbc_param=lst_param, atbc_starting=lst_start))
  
}


atbc.probs <- function(th, param_list, nOccas, sc) {
  
  probs <- rep(NA, nrow(th))
  for (i in 1:nrow(th)) {
    probs[i] <- latent.probs.calc.atbc(nocca=nOccas,
                                       true_hists=th[i,], 
                                       p=param_list[[3]],
                                       phi=param_list[[4]],
                                       gamma=param_list[[5]],
                                       alpha=param_list[[1]],
                                       beta=param_list[[2]],
                                       between.primary.ind=cumsum(sc))
  }
  return(probs)
}

acbt.assign <- function(seed, 
                        prim, 
                        occa, 
                        between.inds, 
                        lst_param=vector(mode="list", length=5), 
                        lst_start=vector(mode="list", length=5),
                        migrationprob,
                        returnprob,
                        captureprob,
                        survivalprob) {
  
  set.seed(seed)
  
  aprim <- migrationprob
  a <- rep(1, occa-1)
  a[cumsum(between.inds)[-prim]] <- aprim
  
  bprim <- rep(returnprob, prim-2)
  b <- rep(1, occa-1)
  b[cumsum(between.inds)[-c(1,prim)]] <- bprim
  
  p <- rep(captureprob, occa)
  
  phiprim <- rep(survivalprob, prim-1)
  phi <- rep(1, occa-1)
  phi[cumsum(between.inds)[-prim]] <- phiprim
  
  gprim <- rep(1/(prim), prim)
  g <- rep(1, occa)
  g[c(1,cumsum(between.inds)+1)[-prim-1]] <- gprim
  
  lst_param[[1]] <- a
  lst_param[[2]] <- b
  lst_param[[3]] <- p
  lst_param[[4]] <- phi
  lst_param[[5]] <- g
  
  lst_start[[1]] <- logit(aprim)
  lst_start[[2]] <- logit(bprim)
  lst_start[[3]] <- logit(p)
  lst_start[[4]] <- logit(phiprim)
  g0 <- gprim[-prim] 
  lst_start[[5]] <- trans.beta.func(g0)
  
  return(list(acbt_param=lst_param, acbt_starting=lst_start))
  
}


acbt.probs <- function(th, param_list, nOccas, sc) {
  
  probs <- rep(NA, nrow(th))
  for (i in 1:nrow(th)) {
    probs[i] <- latent.probs.calc.acbt(nocca=nOccas,
                                       true_hists=th[i,], 
                                       p=param_list[[3]],
                                       phi=param_list[[4]],
                                       gamma=param_list[[5]],
                                       alpha=param_list[[1]],
                                       beta=param_list[[2]],
                                       between.primary.ind=cumsum(sc))
  }
  return(probs)
}

acbc.assign <- function(seed, 
                        prim, 
                        occa, 
                        between.inds, 
                        lst_param=vector(mode="list", length=5), 
                        lst_start=vector(mode="list", length=5),
                        migrationprob,
                        returnprob,
                        captureprob,
                        survivalprob) {
  set.seed(seed)
  
  aprim <- migrationprob
  a <- rep(1, occa-1)
  a[cumsum(between.inds)[-prim]] <- aprim
  
  bprim <- returnprob
  b <- rep(1, occa-1)
  b[cumsum(between.inds)[-c(1,prim)]] <- bprim
  
  p <- rep(captureprob, occa)
  
  phiprim <- rep(survivalprob, prim-1)
  phi <- rep(1, occa-1)
  phi[cumsum(between.inds)[-prim]] <- phiprim
  
  gprim <- rep(1/(prim), prim)
  g <- rep(1, occa)
  g[c(1,cumsum(between.inds)+1)[-prim-1]] <- gprim
  
  lst_param[[1]] <- a
  lst_param[[2]] <- b
  lst_param[[3]] <- p
  lst_param[[4]] <- phi
  lst_param[[5]] <- g
  
  lst_start[[1]] <- logit(aprim)
  lst_start[[2]] <- logit(bprim)
  lst_start[[3]] <- logit(p)
  lst_start[[4]] <- logit(phiprim)
  g0 <- gprim[-prim] 
  lst_start[[5]] <- trans.beta.func(g0)
  
  return(list(acbc_param=lst_param, acbc_starting=lst_start))
  
}

acbc.probs <- function(th, param_list, nOccas, sc) {
  
  probs <- rep(NA, nrow(th))
  for (i in 1:nrow(th)) {
    probs[i] <- latent.probs.calc.acbc(nocca=nOccas,
                                       true_hists=th[i,], 
                                       p=param_list[[3]],
                                       phi=param_list[[4]],
                                       gamma=param_list[[5]],
                                       alpha=param_list[[1]],
                                       beta=param_list[[2]],
                                       between.primary.ind=cumsum(sc))
  }
  return(probs)
}

cr.assign <- function(seed, 
                      prim, 
                      occa, 
                      between.inds, 
                      lst_param=vector(mode="list", length=5), 
                      lst_start=vector(mode="list", length=5),
                      migrationprob,
                      captureprob,
                      survivalprob) {
  
  set.seed(seed)
  
  aprim <- rep(migrationprob, prim-1)
  last.alpha <- aprim[length(aprim)]
  aprim <- aprim[-length(aprim)]
  a <- rep(1, occa-1)
  a[cumsum(between.inds)[1:length(aprim)]] <- aprim
  
  phiprim <- rep(survivalprob, prim-1)
  last.phi <- phiprim[length(phiprim)]
  phiprim <- phiprim[-length(phiprim)]
  phi <- rep(1, occa-1)
  phi[cumsum(between.inds)[-c(prim-1, prim)]] <- phiprim
  
  prod.alpha <- last.phi*(1-last.alpha)
  
  p <- rep(captureprob, occa)
  
  gprim <- rep(1/(prim), prim)
  g <- rep(1, occa)
  g[c(1,cumsum(between.inds)+1)[-prim-1]] <- gprim
  
  lst_param[[1]] <- a
  lst_param[[2]] <- p
  lst_param[[3]] <- phi
  lst_param[[4]] <- g
  lst_param[[5]] <- prod.alpha
  
  lst_start[[1]] <- logit(aprim)
  lst_start[[2]] <- logit(p)
  lst_start[[3]] <- logit(phiprim)
  g0 <- gprim[-prim] 
  lst_start[[4]] <- trans.beta.func(g0)
  lst_start[[5]] <- logit(prod.alpha)
  
  return(list(param=lst_param, starting=lst_start))
  
}


cr.probs <- function(th, param_list, nOccas, sc) {
  
  probs <- rep(NA, nrow(th))
  for (i in 1:nrow(th)) {
    probs[i] <- latent.probs.calc.cr(nocca=nOccas,
                                     true_hists=th[i,], 
                                     p=param_list[[2]],
                                     phi=param_list[[3]],
                                     gamma=param_list[[4]],
                                     alpha=param_list[[1]],
                                     prodalpha=param_list[[5]],
                                     between.primary.ind=cumsum(sc))
  }
  return(probs)
}

