library(RSpectra)
library(MASS)
source("mca.R")
source("paraplu.R")
source("majorization.R")
source("auxiliaries.R")

set.seed(12345)

xx <- matrix(c(
  sample(1:3, 20, replace = TRUE),
  sample(1:5, 20, replace = TRUE),
  sample(1:2, 20, replace = TRUE)
), 20, 3)
gx <- as.list(1:3)
gx[[1]] <- makeIndicator(xx[, 1])
gx[[2]] <- makeIndicator(xx[, 2])
gx[[3]] <- makeIndicator(xx[, 3])

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
  indi <- sapply(1:nvar, list)
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
  dhat <- makeDhat(indi, ncat, norm)
  sold <- sum((dhat - dold)^2)
  repeat {
    haux <- majorizationStep(dold, dhat, xold, yold, vinv, inmax, ips) 
    imtc <- 0L
    sold <- 0.0
    for (i in 1:nobj) {
      ksum <- 0L
      for (j in 1:nvar) {
        kimj <- ncat[j]
        targ <- dold[i, ksum + 1:kimj]
        marg <- mean(targ)
        kind <- which.max(indi[[j]][i, ])
        if (which.min(targ) == kind) {
          imtc <- imtc + 1
          ehat <- targ
        } else {
          ehat <- paraplu(targ, kind)
        }
        dhat[i, ksum + 1:kimj] <- ehat
        ksum <- ksum + kimj
      }
    }
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



