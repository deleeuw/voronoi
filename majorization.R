majorizationStep <- function(dold, dhat, xold, yold, vinv, inmax, ips) {
  bmat <- makeBmat(dold, dhat)
  xmid <- bmat$b11 * xold + bmat$b12 %*% yold
  ymid <- bmat$b22 * yold + crossprod(bmat$b12, xold)
}

makeBmat <- function(dold, dhat) {
  ktot <- ncol(dold)
  nobj <- nrow(dold)
  b12 <- -dhat / dold
  b11 <- -rowSums(b12)
  b22 <- -colSums(b12)
  return(list(b11 = b11, b22 = b22, b12 = b12))
}

makeVinv <- function(n, m) {
  vinv11 <- ((diag(n) - (1 / n)) / m) + ((m / n) / ((n + m)^2))
  vinv12 <- matrix(-1, n, m) / ((n + m)^2)
  vinv22 <- ((diag(m) - (1 / m)) / n) + ((n / m) / ((n + m)^2))
  return(list(vinv11 = vinv11, vinv22 = vinv22, vinv12 = vinv12))
}