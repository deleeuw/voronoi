mca <- function(indi, ndim) {
  marg <- colSums(indi)
  pmat <- makePmat(indi)
  eigs <- eigs_sym(pmat, ndim + 1)
  xmat <- eigs$vectors[, -1]
  ymat <- crossprod(indi, xmat) / marg
  return(list(xmat = xmat, ymat = ymat))
}

makePmat <- function(indi) {
  marg <- 1 / colSums(indi)
  return(tcrossprod(indi %*% diag(marg), indi))
}