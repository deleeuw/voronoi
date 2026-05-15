
monotone <- function(dold, ncat, indi) {
  nobj <- nrow(dold)
  nvar <- length(ncat)
  dhat <- array(0, dim(dold))
  for (i in 1:nobj) {
    ksum <- 0L
    for (j in 1:nvar) {
      kimj <- ncat[j]
      targ <- dold[i, ksum + 1:kimj]
      marg <- mean(targ)
      kind <- which.max(indi[i, ksum + 1:kimj])
      if (which.min(targ) == kind) {
        ehat <- targ
      } else {
        ehat <- paraplu(targ, kind)
      }
      dhat[i, ksum + 1:kimj] <- ehat
      ksum <- ksum + kimj
    }
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

normalizeDhat <- function(dhat, dnew, indi, ncat, dnorm) {
  nvar <- length(ncat)
  nobj <- nrow(dhat)
  kwrk <- 0
  for (j in 1:nvar) {
    for (i in 1:nobj) {
      nwrk <- ncat[j]
      hwrk <- dhat[i, kwrk + 1:nwrk]
      dave <- mean(hwrk)
      dmax <- max(abs(hwrk - dave))
      if (dmax < 1e-10) {
        iwrk <- which(indi[i, kwrk + 1:nwrk] == 1)
        hwrk[iwrk] <- -Inf
        jwrk <- which.max(hwrk)
        hwrk <- ifelse(jwrk == 1:nwrk, nwrk, 0) - 1
        dhat[i, kwrk + 1:nwrk] <- hwrk / sqrt(nwrk * (nwrk - 1)) + dave
      } else {
        hwrk <- hwrk - dave
        dhat[i, kwrk + 1:nwrk] <- hwrk / sqrt(sum(hwrk^2)) + dave
      }
    }
    kwrk <- kwrk + ncat[j]
  }
  return(dhat)
}