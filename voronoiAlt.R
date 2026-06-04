voronoiHomogeneityAnalysis <- function(x,
                                       ndim = 2,
                                       wght = NULL,
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
  indi <- rep(list(0), nvar)
  for (j in 1:nvar) {
    indi[[j]] <- makeIndicator(x[, j])
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
  xold <- haux$xmat
  yold <- haux$ymat
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
    bmat <- rep(list(0), nvar)
    for (j in 1:nvar) {
      bmat[[j]] <- makeBmat(wght[[j]], dold[[j]], dhat[[j]])
    }
    bsum <- rowSums(sapply(bmat, rowSums))
    xnew <- bsum * xold
    for (j in 1:nvar) {
      xnew <- xnew - (bmat[[j]] - wght[[j]]) %*% yold[[j]]
    }
    xnew <- xnew / wsum
    dnew <- rep(list(0), nvar)
    ynew <- rep(list(0), nvar)
    snew <- 0
    daps <- 0
    for (j in 1:nvar) {
      dnew[[j]] <- makeDmat(xnew, yold[[j]])
      bmat[[j]] <- makeBmat(wght[[j]], dnew[[j]], dhat[[j]])
      ynew[[j]] <- colSums(bmat[[j]]) * yold[[j]]
      ynew[[j]] <- ynew[[j]] - crossprod(bmat[[j]] - wght[[j]], xnew)
      ynew[[j]] <- ynew[[j]] / colSums(wght[[j]])
      dnew[[j]] <- makeDmat(xnew, ynew[[j]])
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
  return(list(
    xmat = xnew,
    ymat = ynew,
    dmat = dnew,
    dhat = dhat,
    itel = itel,
    loss = snew,
    indi = indi,
    wght = wght
  ))
}

mca <- function(indi,
                ncat,
                ndim,
                itmax = 1000,
                eps = 1e-6,
                verbose = FALSE) {
  nobj <- nrow(indi[[1]])
  nvar <- length(indi)
  xold <- matrix(rnorm(nobj * ndim), nobj, ndim)
  xold <- apply(xold, 2, function(x)
    x - mean(x))
  xold <- qr.Q(qr(xold))
  ymat <- rep(list(0), nvar)
  itel <- 1
  repeat {
    xnew <- matrix(0, nobj, ndim)
    for (j in 1:nvar) {
      ymat[[j]] <- crossprod(indi[[j]], xold) / colSums(indi[[j]])
      xnew <- xnew + indi[[j]] %*% ymat[[j]]
    }
    xnew <- apply(xnew, 2, function(x)
      x - mean(x))
    xnew <- qr.Q(qr(xnew))
    aps <- max(abs(xold - xnew))
    if (verbose) {
      cat(
        "itel",
        formatC(itel, digits = 4, format = "d"),
        "aps ",
        formatC(
          aps,
          digits = 10,
          width = 15,
          format = "f"
        ),
        "\n"
      )
    }
    if ((itel == itmax) || (aps < eps)) {
      break
    }
    xold <- xnew
    itel <- itel + 1
  }
  return(list(xmat = xnew, ymat = ymat))
}

makeIndicator <- function(x, labs = NULL) {
  return(ifelse(outer(x, sort(unique(
    x
  )), "=="), 1, 0))
}

makeDmat <- function(x, y) {
  xx <- rowSums(x^2)
  yy <- rowSums(y^2)
  dmat <- sqrt(outer(xx, yy, "+") - 2 * tcrossprod(x, y))
  return(dmat)
}

makeBmat <- function(wght, dmat, dhat) {
  return(wght * dhat / dmat)
}

monotone <- function(dold, ncat, indi) {
  nobj <- nrow(dold)
  dhat <- matrix(0, nobj, ncat)
  for (i in 1:nobj) {
    targ <- dold[i, ]
    kind <- which.max(indi[i, ])
    if (which.min(targ) == kind) {
      ehat <- targ
    } else {
      ehat <- paraplu(targ, kind)
    }
    dhat[i, ] <- ehat
  }
  return(dhat)
}

paraplu <- function(x, k) {
  n <- length(x)
  tx <- rep(0, n)
  xk <- x[k]
  if (xk == min(x)) {
    return(x)
  }
  x[k] <- -1
  h <- sort(x, index.return = TRUE)
  sx <- h$x
  sx[1] <- xk
  ix <- h$ix
  for (i in 2:n) {
    r <- mean(sx[1:i])
    if (i == n) {
      return(rep(r, n))
    }
    if (r < sx[i + 1]) {
      tx <- c(rep(r, i), sx[(i + 1):n])
      break
    }
  }
  for (i in 1:n) {
    sx[ix[i]] <- tx[i]
  }
  return(sx)
}

vha <- voronoiHomogeneityAnalysis
