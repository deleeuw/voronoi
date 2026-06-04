ei <- function(i, n) {
  ifelse(i == 1:n, 1, 0)
}

aij <- function(i, j, n) {
  e1 <- ei(i, n)
  e2 <- ei(j, n)
  return(outer(e1 - e2, e1 - e2))
}

makeMatrixFromList <- function(y) {
  ktot <- sum(sapply(y, ncol))
  nobj <- nrow(y[[1]])
  return(matrix(unlist(y), nobj, ktot))
}

makeListfromMatrix <- function(y, k) {
  nmat <- length(k)
  nobj <- nrow(y)
  ylist <- list(1:nmat)
  ksum <- 0
  for (j in 1:nmat) {
    ylist[[j]] <- y[ , ksum + 1:k[j]]
    ksum <- ksum + k[j]
  }
  return(ylist)
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

