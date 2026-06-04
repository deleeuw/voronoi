suppressPackageStartupMessages(library(RSpectra, quietly = TRUE))
suppressPackageStartupMessages(library(deldir, quietly = TRUE))
suppressPackageStartupMessages(library(dismo, quietly = TRUE))

source("mca.R")
source("majorization.R")
source("auxiliaries.R")
source("monotone.R")

voronoiHomogeneityAnalysis <- function(x,
                                       ndim = 2,
                                       inmax = 5,
                                       ips = 1e-6,
                                       itmax = 1000,
                                       eps = 1e-6,
                                       aps = 1e-6,
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
  if (length(yrank) == 1) {
    yrank <- rep(yrank, nvar)
  }
  ktot <- ncol(indi)
  vinv <- makeVinv(nobj, ktot)
  hini <- mca(indi, ndim)
  xold <- hini$xmat
  yold <- hini$ymat
  dold <- makeDmat(xold, yold)
  dhat <- monotone(dold, ncat, indi)
  if (dnorm) {
    dhat <- normalizeDhat(dhat, dold, indi, ncat, dnorm)
  }
  sold <- sum((dhat - dold)^2)
  itel <- 1
  repeat {
    haux <- majorizationStep(dold, dhat, xold, yold, vinv, inmax, ips)
    xnew <- haux$xnew
    ynew <- haux$ynew
    meax <- apply(xnew, 2, mean)
    xnew <- xnew - outer(rep(1, nobj), meax)
    # ynew <- ynew - outer(rep(1, ktot), meax)
    if (xnorm) {
      xsvd <- svd(xnew)
      xnew <- tcrossprod(xsvd$u, xsvd$v)
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
    if (dnorm) {
      dhat <- normalizeDhat(dhat, dnew, indi, ncat, dnorm)
    }
    snew <- sum((dhat - dnew)^2)
    daps <- max(abs(dold - dnew))
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
        "daps",
        formatC(
          daps,
          digits = 10,
          width = 15,
          format = "f"
        ),
        "\n"
      )
    }
    if ((itel == itmax) || (daps < eps) || (snew < aps)) {
      break
    }
    dold <- dnew
    xold <- xnew
    yold <- ynew
    sold <- snew
    itel <- itel + 1
  }
  return(
    list(
      xmat = xnew,
      ymat = ynew,
      dmat = dnew,
      dhat = dhat,
      itel = itel,
      loss = snew,
      xini = hini$xmat,
      yini = hini$ymat,
      indi = indi,
      ncat = ncat
    )
  )
}

vha <- voronoiHomogeneityAnalysis

