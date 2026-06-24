library(dismo)

plotVHA <- function(xmat,
                    ymat,
                    indi,
                    main = NULL,
                    ylab = NULL,
                    xlab = 1,
                    pdf = FALSE) {
  ncat <- ncol(indi)
  rmin <- 1.1 * min(rbind(xmat, ymat))
  rmax <- 1.1 * max(rbind(xmat, ymat))
  if (pdf) {
    pdf(paste("plot", main, ".pdf", sep = ""))
  }
  plot(voronoi(ymat, ext = c(rmin, rmax, rmin, rmax)), main = main)
  if (length(xlab) > 1) {
    text(xmat, as.character(xlab), col = "BLACK")
  } else {
    if (xlab == 0) {
      text(xmat, as.character(1:nrow(xmat)), col = "BLACK")
    }
    if (xlab == 1) {
      text(xmat, as.character(drop(indi %*% 1:ncat)), col = "BLACK")
    }
  }
  if (!is.null(ylab)) {
    text(ymat, ylab, col = "BLACK", cex = 1.5)
  } else {
    points(ymat,
           col = "RED",
           cex = 1.5,
           pch = 16)
  }
  if (pdf) {
    dev.off()
  }
}
