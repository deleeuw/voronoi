
binaryMonotoneRegression <- function(y, w, g, verbose = TRUE) {
  nn <- length(y)
  x <- rep(0, nn)
  n0 <- which(g == 0)
  n1 <- which(g == 1)
  sy <- sort(y)
  gold <- -Inf
  for (k in 1:nn) {
    anew <- sy[k]
    x0 <- pmin(anew - y[n0], 0)
    x1 <- pmax(anew - y[n1], 0)
    gnew <- sum(w[n0] * x0) + sum(w[n1] * x1)
    fnew <- sum(w[n0] * x0^2) + sum(w[n1] * x1^2)
    if (verbose) {
      cat("alph =", formatC(anew, digits = 6, width = 10, format = "f"),
          "loss =", formatC(fnew, digits = 6, width = 10, format = "f"),
          "grad =", formatC(gnew, digits = 6, width = 10, format = "f"),
          "\n")
    }
    if (gnew >= 0) {
      break
    }
    gold <- gnew
    aold <- anew
  }
  if (gnew == 0) {
    ao <- anew
  } else {
    as <- (gnew - gold) / (anew - aold)
    ao <- anew - gnew / as
  }
  x[n0] <- pmin(y[n0], ao)
  x[n1] <- pmax(y[n1], ao)
  return(x)
}

bmr <- binaryMonotoneRegression