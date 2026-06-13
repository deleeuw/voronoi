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
                                       xnorm = FALSE,
                                       yrank = ndim,
                                       verbose = TRUE) {
  nobj <- nrow(theData)
  nvar <- ncol(theData)
  ncat <- rep(0, nvar)
  indi <- rep(list(0), nvar)
  for (j in 1:nvar) {
    indi[[j]] <- makeIndicator(theData[, j])
    ncat[j] <- ncol(indi[[j]])
  }
  if (is.null(wght)) {
    wght <- rep(list(0), nvar)
    for (j in 1:nvar) {
      wght[[j]] <- matrix(1, nobj, ncat[j])
    }
  }
  wsum <- rowSums(sapply(wght, rowSums))
  haux <- mca(indi, ncat, ndim)
  xini <- xold <- haux$xmat
  yini <- yold <- haux$ymat
  dold <- rep(list(0), nvar)
  dhat <- rep(list(0), nvar)
  sold <- 0.0
  for (j in 1:nvar) {
    dold[[j]] <- makeDmat(xold, yold[[j]])
    dhat[[j]] <- monotone(dold[[j]], ncat[j], indi[[j]])
    sold <- sold + sum(wght[[j]] * (dhat[[j]] - dold[[j]])^2)
  }
  itel <- 1
  repeat {
    xnew <- makeXnew(xold, yold, dold, dhat, wght, wsum, xcent, xnorm)
    ynew <- makeYnew(xnew, yold, dold, dhat, wght, yrank)
    snew <- 0.0
    smid <- 0.0
    daps <- 0.0
    dnew <- rep(list(0), nvar)
    for (j in 1:nvar) {
      dnew[[j]] <- makeDmat(xnew, ynew[[j]])
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
