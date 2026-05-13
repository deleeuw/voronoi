library(RSpectra)
library(deldir)
library(dismo)

source("mca.R")
source("paraplu.R")
source("majorization.R")
source("auxiliaries.R")
source("monotone.R")

voronoiHomogeneityAnalysis <- function(x,
                                       ndim = 2,
                                       inmax = 5,
                                       ips = 1e-6,
                                       itmax = 1000,
                                       eps = 1e-6,
                                       dnorm = FALSE,
                                       xnorm = FALSE,
                                       yrank = NULL,
                                       verbose = FALSE) {
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
  xold <- haux$xmat
  yold <- haux$ymat
  dold <- makeDmat(xold, yold)
  # initialize dhat from indi
  # dhat <- initializeDhat
  dhat <- monotone(dold, ncat, indi)
  sold <- sum((dhat - dold)^2)
  itel <- 1
  repeat {
    haux <- majorizationStep(dold, dhat, xold, yold, vinv, inmax, ips)
    xnew <- haux$xnew
    ynew <- haux$ynew
    meax <- apply(xnew, 2, mean)
    xnew <- xnew - outer(rep(1, nobj), meax)
    ynew <- ynew - outer(rep(1, ktot), meax)
    if (xnorm) {
      xnew <- qr.Q(qr(xnew))
    }
    if (!is.null(yrank)) {
      ksum <- 0
      for (j in 1:nvar) {
        kind <- ksum + 1:ncat[j]
        if (yrank[j]) {
          hsvd <- svd(ynew[kind, ], nu = 1, nv = 1)
          ynew[kind, ] <- hsvd$d[1] * outer(drop(hsvd$u), drop(hsvd$v))
        }
        ksum <- ksum + ncat[j]
      }
    }
    dnew <- makeDmat(xnew, ynew)
    smid <- sum((dhat - dnew)^2)
    dhat <- monotone(dnew, ncat, indi)
    snew <- sum((dhat - dnew)^2)
    if (verbose) {
      cat(
        "itel",
        formatC(itel, digits = 4, format = "d"),
        "sold",
        formatC(
          sold,
          digits = 10,
          width = 15,
          format = "f"
        ),
        "smid",
        formatC(
          smid,
          digits = 10,
          width = 15,
          format = "f"
        ),
        "snew",
        formatC(
          snew,
          digits = 10,
          width = 15,
          format = "f"
        ),
        "\n"
      )
    }
    aps <- max(abs(dold - dnew))
    if ((itel == itmax) || (aps < eps)) {
      break
    }
    dold <- dnew
    xold <- xnew
    yold <- ynew
    sold <- snew
    itel <- itel + 1
  }
  return(list(
    xmat = xnew,
    ymat = ynew,
    dmat = dnew,
    dhat = dhat,
    itel = itel,
    loss = snew,
    xini = haux$xmat,
    yini = haux$ymat,
    indi = indi,
    ncat = ncat
  ))
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

vha <- voronoiHomogeneityAnalysis
