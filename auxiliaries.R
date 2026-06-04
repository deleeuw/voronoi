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