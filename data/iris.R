data(iris)

imat <- matrix(0, 150, 4)
for (j in 1:4) {
  q <- quantile(iris[,j], seq(.05, 1, .25))
  g <- ifelse(outer(iris[, j], q, "<"), 1, 0)
  imat[, j] <- 1 + t(apply(g, 1, diff)) %*% 1:(length(q) - 1)
}

# h<-vha(cbind(imat, iris$Species), verbose = TRUE, yrank = c(1,1,1,1,2))
# plotVHA(h$xini, h$yini[[5]], h$indi[[5]], h$ncat[5], xlab = as.character(iris[, 5]))
# plotVHA(h$xini, h$yini[[4]], h$indi[[4]], h$ncat[4], ylab = as.character(1:4), xlab = as.character(iris[, 5]))
# plotVHA(h$xmat, h$ymat[[5]], h$indi[[5]], h$ncat[5], xlab = as.character(iris[, 5]))
# plotVHA(h$xmat, h$ymat[[4]][1:2, ], h$indi[[4]], ncat = 2, ylab = c("1 4", "2 3"), xlab = as.character(iris[, 5]))