


library(dismo)
library(plotrix)
library(grDevices)

plotVHA <- function(xmat,
                    ymat,
                    indi,
                    dmat,
                    type = "voronoi",
                    main = NULL,
                    ylab = 1,
                    xlab = 1,
                    pdf = FALSE) {
  ncat <- ncol(indi)
  rmin <- 1.1 * min(rbind(xmat, ymat))
  rmax <- 1.1 * max(rbind(xmat, ymat))
  if (pdf) {
    pdf(paste("plot", main, ".pdf", sep = ""))
  }
  plot(
    0,
    xlim = c(rmin, rmax),
    ylim = c(rmin, rmax),
    type = "n",
    main = main
  )
  if (type == "voronoi") {
    plot(voronoi(ymat, ext = c(rmin, rmax, rmin, rmax)), border = "BLUE", lwd = 2)
  }
  if (type == "hull") {
    for (j in 1:ncat) {
      xhul <- xmat[which(indi[, j] == 1), ]
      polygon(xhul[chull(xhul), ], border = "BLUE", lwd = 2)
    }
  }
  if (type == "circle") {
    for (j in 1:ncat) {
      r <- max(dmat[which(indi[, j] == 1), j])
      draw.circle(ymat[j, 1], ymat[j, 2], r, border = "BLUE", lwd = 2)
    }
  }
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
  if (length(ylab) > 1) {
    text(ymat, ylab, col = "BLACK", cex = 1.5)
  } else {
    if (ylab == 0) {
      points(ymat,
             col = "RED",
             cex = 1.5,
             pch = 16)
    }
    if (ylab == 1) {
      text(ymat,
           as.character(1:ncat),
           col = "RED",
           cex = 1.5)
    }
  }
  if (pdf) {
    dev.off()
  }
}
