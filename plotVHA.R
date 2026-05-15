library(dismo)

plotVHA <- function(xmat,
                    ymat,
                    ncat,
                    main = NULL,
                    labs = NULL,
                    pdf = FALSE) {
  nvar <- length(ncat)
  kcat <- 0
  rmin <- 1.1 * min(rbind(xmat, ymat))
  rmax <- 1.1 * max(rbind(xmat, ymat))
  for (j in 1:nvar) {
    if (pdf) {
      pdf(paste("plot", j , ".pdf", sep = ""))
    }
    y <- ymat[kcat + 1:ncat[j], ]
    plot(voronoi(y, ext = c(rmin, rmax, rmin, rmax)), main = main[j])
    text(xmat, as.character(1:nrow(xmat)))
    if (!is.null(labs)) {
      ylab <- labs[kcat + 1:ncat[j]]
      text(y, ylab, col = "RED", cex = 1.5)
    } else {
      points(y,
             col = "RED",
             cex = 1.5,
             pch = 16)
    }
    kcat <- kcat + ncat[j]
    if (pdf) {
      dev.off()
    }
  }
}