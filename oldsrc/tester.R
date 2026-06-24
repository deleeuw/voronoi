n <- 20
i <- 10
l <- 4
p <- 2
#set.seed(12345)
w <- matrix(rnorm(80)^2, n, l)
rw <- rowSums(w)
cw <- colSums(w)
a <- sample(1:l, n, replace = TRUE)
g <- ifelse(outer(a, 1:l, "=="), 1, 0)
d <- colSums(g)
x <- matrix(rnorm(n * p), n, p)
y <- crossprod(g, x) / d
e <- outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y)
# f <- ifelse(i == 1:n, 1, 0)
# h <- outer(f, f) - (outer(g[, l], f) + outer(f, g[, l])) / d[l]
# h <- h + outer(g[, l], g[, l]) / d[l]^2
# v <- sum(x * (h %*% x))
# r <- f - g[, l] / d[l]
# s <- sum(colSums(r * x)^2)
# print(c(e[i, l], v, s))
k <- diag(rw) - g %*% (t(w) / d) - w %*% (t(g) / d) + g %*% (t(g) * (cw / d^2))
print(c(sum(w * e), sum(x * (k %*% x))))

