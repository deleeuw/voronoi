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

makeYnew <- function(xnew, yold, dold, dhat, wght, yrank) {
  nvar <- length(yold)
  ndim <- ncol(xnew)
  dnew <- rep(list(0), nvar)
  ynew <- rep(list(0), nvar)
  bmat <- rep(list(0), nvar)
  if (length(yrank) == 1) {
    yrank <- rep(yrank, nvar)
  }
  wsum <- 
  snew <- 0
  daps <- 0
  for (j in 1:nvar) {
    wsum <- colSums(wght[[j]])
    yrnk <- yrank[j]
    dnew[[j]] <- makeDmat(xnew, yold[[j]])
    bmat[[j]] <- makeBmat(wght[[j]], dnew[[j]], dhat[[j]])
    ynew[[j]] <- colSums(bmat[[j]]) * yold[[j]]
    ynew[[j]] <- ynew[[j]] - crossprod(bmat[[j]] - wght[[j]], xnew)
    ynew[[j]] <- ynew[[j]] / wsum
    if (yrnk < ndim) {
      saux <- svd(sqrt(wsum) * ynew[[j]], nu = yrnk, nv = yrnk)
      daux <- saux$d[1:yrnk]
      daux <- ifelse(yrnk == 1, daux, diag(daux))
      ynew[[j]] <- tcrossprod((saux$u %*% daux) / sqrt(wsum), saux$v)
    }
  }
  return(ynew)
}
