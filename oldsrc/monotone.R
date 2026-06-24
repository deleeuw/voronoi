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
