#===============================================================================
# helpers.R  -  small functions used by all other scripts
#===============================================================================

logit <- binomial()$linkfun

ilogit <- function(x) (1 + exp(-x))^-1

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

catline <- function(...) cat(..., "\n")

# First and last secondary occasion of every primary period
find.primperiod.indices <- function(nOccas) {
  last  <- cumsum(nOccas)
  first <- last - nOccas + 1
  ind <- cbind(first, last)
  dimnames(ind) <- list(paste("Period", seq_along(nOccas)),
                        c("Starting occasion", "Ending occasion"))
  ind
}

# Wald critical value used for every confidence interval (as in the original code)
Z.CRIT <- 1.959964
