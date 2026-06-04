makeXnew <- function(xold, yold, dold, dhat, wght, wsum, xcent, xnorm) {
  nobj <- nrow(xold)
  ndim <- ncol(xold)
  nvar <- length(yold)
  bmat <- rep(list(0), nvar)
  xnew <- matrix(0, nobj, ndim)
  for (j in 1:nvar) {
    bmat[[j]] <- makeBmat(wght[[j]], dold[[j]], dhat[[j]])
    xnew <- xnew + (bmat[[j]] - wght[[j]]) %*% yold[[j]]
  }
  bsum <- rowSums(sapply(bmat, rowSums))
  xnew <- bsum * xold - xnew
  xnew <- xnew / wsum
  if (xcent) {
    xnew <- xnew - outer(rep(1, nobj), drop(wsum %*% xnew)) / sum(wsum)
  }
  if (xnorm) {
    xaux <- sqrt(wsum) * xnew
    saux <- svd(xaux)
    xnew <- tcrossprod(saux$u , saux$v) / sqrt(wsum)
  }
  return(xnew)
}

makeYnew <- function(xnew, yold, dold, dhat, wght) {
  nvar <- length(yold)
  dnew <- rep(list(0), nvar)
  ynew <- rep(list(0), nvar)
  bmat <- rep(list(0), nvar)
  snew <- 0
  daps <- 0
  for (j in 1:nvar) {
    dnew[[j]] <- makeDmat(xnew, yold[[j]])
    bmat[[j]] <- makeBmat(wght[[j]], dnew[[j]], dhat[[j]])
    ynew[[j]] <- colSums(bmat[[j]]) * yold[[j]]
    ynew[[j]] <- ynew[[j]] - crossprod(bmat[[j]] - wght[[j]], xnew)
    ynew[[j]] <- ynew[[j]] / colSums(wght[[j]])
  }
  return(ynew)
}
