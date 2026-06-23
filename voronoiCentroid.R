source("mca.R")
source("mals.R")
source("monotone.R")
source("auxiliaries.R")

voronoiCentroidAnalysis <- function(theData,
                                    ndim = 2,
                                    wght = NULL,
                                    inmax = 5,
                                    ips = 1e-6,
                                    itmax = 1000,
                                    eps = 1e-6,
                                    aps = 1e-6,
                                    dnorm = FALSE,
                                    xcent = FALSE,
                                    xnorm = 1,
                                    yrank = ndim,
                                    verbose = TRUE) {
  nobj <- nrow(theData)
  nvar <- ncol(theData)
  indi <- lapply(theData, makeIndicator)
  marg <- lapply(indi, colSums)
  ncat <- sapply(indi, ncol)
  if (is.null(wght)) {
    wght <- lapply(indi, function(x)
      array(1, dim(x)))
  }
  wrsm <- lapply(wght, rowSums)
  wcsm <- lapply(wght, colSums)
  vmat <- matrix(0, nobj, nobj)
  for (j in 1:nvar) {
    vmat <- vmat + diag(wrsm[[j]])
    vmat <- vmat - indi[[j]] %*% (t(wght[[j]]) / marg[[j]])
    vmat <- vmat - wght[[j]] %*% (t(indi[[j]]) / marg[[j]])
    vmat <- vmat + indi[[j]] %*% (t(indi[[j]]) * (wcsm[[j]] / marg[[j]]^2))
  }
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  haux <- mca(indi, ncat, ndim)
  xini <- xold <- haux$xmat
  yini <- lapply(indi, function(x)
    crossprod(x, xini) / colSums(x))
  yold <- yini
  dhat <- rep(list(0), nvar)
  dold <- lapply(yold, function(x)
    makeDmat(xold, x))
  sold <- 0.0
  for (j in 1:nvar) {
    dhat[[j]] <- monotone(dold[[j]], ncat[j], indi[[j]])
    sold <- sold + sum(wght[[j]] * (dhat[[j]] - dold[[j]])^2)
  }
  itel <- 1
  repeat {
    bmat <- matrix(0, nobj, nobj)
    for (j in 1:nvar) {
      resi <- wght[[j]] * dhat[[j]] * invMe(dold[[j]])
      rrsm <- rowSums(resi)
      rcsm <- colSums(resi)
      bmat <- bmat + diag(rrsm)
      bmat <- bmat - indi[[j]] %*% (t(resi) / marg[[j]])
      bmat <- bmat - resi %*% (t(indi[[j]]) / marg[[j]])
      bmat <- bmat + indi[[j]] %*% (t(indi[[j]]) * (rcsm / marg[[j]]^2))
    }
    xnew <- vinv %*% bmat %*% xold
    xnew <- switch(xnorm, 
                   xnew,
                   xnew %*% matrixPower(crossprod(xnew, vmat %*% xnew), -0.5),
                   xnew / sqrt(sum(xnew * (vmat %*% xnew))))
    ynew <- lapply(indi, function(x)
     crossprod(x, xnew) / colSums(x))
    snew <- 0.0
    smid <- 0.0
    daps <- 0.0
    dnew <- lapply(ynew, function(x)
      makeDmat(xnew, x))
    for (j in 1:nvar) {
      smid <- smid + sum(wght[[j]] * (dhat[[j]] - dnew[[j]])^2)
      dhat[[j]] <- monotone(dnew[[j]], ncat[j], indi[[j]])
      snew <- snew + sum(wght[[j]] * (dhat[[j]] - dnew[[j]])^2)
      daps <- max(daps, max(abs(dold[[j]] - dnew[[j]])))
    }
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
    if ((itel == itmax) || (daps < eps)) {
      break
    }
    xold <- xnew
    yold <- ynew
    dold <- dnew
    sold <- snew
    itel <- itel + 1
  }
  return(
    list(
      xmat = xnew,
      ymat = ynew,
      dmat = dnew,
      dhat = dhat,
      xini = xini,
      yini = yini,
      itel = itel,
      loss = snew,
      indi = indi,
      ncat = ncat,
      wght = wght
    )
  )
}

vca <- voronoiCentroidAnalysis
