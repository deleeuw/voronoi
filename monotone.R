source("paraplu.R")

monotone <- function(dhat, dold, ncat, indi) {
  nobj <- nrow(dhat)
  nvar <- length(ncat)
  imtc <- 0L
  for (i in 1:nobj) {
    ksum <- 0L
    for (j in 1:nvar) {
      kimj <- ncat[j]
      targ <- dold[i, ksum + 1:kimj]
      marg <- mean(targ)
      kind <- which.max(indi[i, ksum + 1:kimj])
      if (which.min(targ) == kind) {
        imtc <- imtc + 1
        ehat <- targ
      } else {
        ehat <- paraplu(targ, kind)
      }
      dhat[i, ksum + 1:kimj] <- ehat
      ksum <- ksum + kimj
    }
  }
}
