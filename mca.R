mca <- function(indi, ndim) {
  nobj <- nrow(indi[[1]])
  marg <- colSums(indi)
  pmat <- makePmat(indi)
  eigs <- eigs_sym(pmat, ndim + 1)
  xold <- eigs$vectors[, -1]
  yold <- crossprod(indi, xold) / marg
  return(list(xold = xold, yold = yold))
}

makePmat <- function(indi) {
  marg <- 1 / colSums(indi)
  return(tcrossprod(indi %*% diag(marg), indi))
}