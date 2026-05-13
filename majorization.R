majorizationStep <- function(dold, dhat, xold, yold, vinv, inmax, ips) {
  n <- nrow(xold)
  m <- nrow(yold)
  bmat <- makeBmat(dold, dhat)
  xmid <- bmat$b11 * xold + bmat$b12 %*% yold
  ymid <- bmat$b22 * yold + crossprod(bmat$b12, xold)
  vinv <- makeVinv(n, m)
  xnew <- vinv$v11 %*% xmid + vinv$v12 %*% ymid
  ynew <- t(vinv$v12) %*% xmid + vinv$v22 %*% ymid
  return(list(xnew = xnew, ynew = ynew))
}

makeBmat <- function(dold, dhat) {
  ktot <- ncol(dold)
  nobj <- nrow(dold)
  b12 <- -dhat / ifelse(dold < 1e-10, 1, dold)
  b11 <- -rowSums(b12)
  b22 <- -colSums(b12)
  return(list(b11 = b11, b22 = b22, b12 = b12))
}

makeVinv <- function(n, m) {
  v11 <- ((diag(n) - (1 / n)) / m) + ((m / n) / ((n + m)^2))
  v12 <- matrix(-1, n, m) / ((n + m)^2)
  v22 <- ((diag(m) - (1 / m)) / n) + ((n / m) / ((n + m)^2))
  return(list(v11 = v11, v22 = v22, v12 = v12))
}