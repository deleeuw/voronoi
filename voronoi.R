library(RSpectra)
library(MASS)
source("mca.R")
source("paraplu.R")
source("majorization.R")
source("auxiliaries.R")


mdshom <- function(x,
                   ndim = 2,
                   norm = 1,
                   inmax = 5,
                   ips = 1e-6,
                   itmax = 100,
                   eps = 1e-6,
                   verbose = TRUE) {
  nobj <- nrow(x)
  nvar <- ncol(x)
  ncat <- rep(0, nvar)
  indi <- NULL
  for (j in 1:nvar) {
    indj <- makeIndicator(x[, j])
    indi <- cbind(indi, indj)
    ncat[j] <- ncol(indj)
  }
  ktot <- ncol(indi)
  vinv <- makeVinv(nobj, ktot)
  haux <- mca(indi, ndim)
  xold <- haux$xold
  yold <- haux$yold
  dold <- makeDmat(xold, yold)
  dhat <- monotone(dhat, dold, ncat, indi)
  sold <- sum((dhat - dold)^2)
  repeat {
    haux <- majorizationStep(dold, dhat, xold, yold, vinv, inmax, ips) 
    imtc <- 0L
    sold <- 0.0
    sold <- sum((dhat - dold)^2)
    mold <- imtc / (nvar * nobj)
    print(c(mold, sold))
    print(dhat)
    break
  }
}

makeDmat <- function(x, y) {
  xx <- rowSums(x^2)
  yy <- rowSums(y^2)
  dmat <- sqrt(outer(xx, yy, "+") - 2 * tcrossprod(x, y))
  return(dmat)
}

makeDhat <- function(indi, ncat, norm = 1) {
  nobj <- nrow(indi)
  nvar <- length(ncat)
  dhat <- NULL
  ksum <- 0
  for (j in 1:nvar) {
    catj <- ncat[j]
    fact <- switch(norm, sqrt(catj / (catj - 1)), sqrt(catj^2 / (catj - 1)))
    indj <- indi[, ksum + 1:catj]
    dhat <- cbind(dhat, indj * fact)
  }
  return(dhat)
}



