
mca <- function(indi,
                ncat,
                ndim,
                itmax = 1000,
                eps = 1e-6,
                verbose = FALSE) {
  nobj <- nrow(indi[[1]])
  nvar <- length(indi)
  xold <- matrix(rnorm(nobj * ndim), nobj, ndim)
  xold <- apply(xold, 2, function(x)
    x - mean(x))
  xold <- qr.Q(qr(xold))
  ymat <- rep(list(0), nvar)
  itel <- 1
  repeat {
    xnew <- matrix(0, nobj, ndim)
    for (j in 1:nvar) {
      ymat[[j]] <- crossprod(indi[[j]], xold) / colSums(indi[[j]])
      xnew <- xnew + indi[[j]] %*% ymat[[j]]
    }
    xnew <- apply(xnew, 2, function(x)
      x - mean(x))
    xnew <- qr.Q(qr(xnew))
    aps <- max(abs(xold - xnew))
    if (verbose) {
      cat(
        "itel",
        formatC(itel, digits = 4, format = "d"),
        "aps ",
        formatC(
          aps,
          digits = 10,
          width = 15,
          format = "f"
        ),
        "\n"
      )
    }
    if ((itel == itmax) || (aps < eps)) {
      break
    }
    xold <- xnew
    itel <- itel + 1
  }
  return(list(xmat = xnew, ymat = ymat))
}
