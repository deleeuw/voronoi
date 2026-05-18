library(RSpectra)
indep 
indep <- function(a, k) {
  x <- drop(eigs_sym(a, 1)$vectors)
  n <- nrow(a)
  q <- c(quantile(x, probs = (0:k) / k), Inf)
  g <- ifelse(outer(x, q, "<"), 1, 0)
  h <- matrix(0, n, k + 1)
  for (j in 1:(k + 1)) {
    h[, j] <- g[, j + 1] - g[, j]
  }
  return(h)
}
  