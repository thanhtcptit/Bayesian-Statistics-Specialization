rm(list = ls())
set.seed(81196)

dat = read.csv("../dataset/nestsize.csv", header = FALSE)
x = dat$V1
n = length(x)

KK = 2
w = 1 / 2
lambda = mean(x)

xx = seq(0, max(x))
yy = w * (xx == 0) + (1 - w) * dpois(xx, lambda)
par(mfrow = c(2, 1))
plot(table(x), type = "h")
plot(xx, yy, type = "h", ylim = c(0, max(yy)), xlab = "x", ylab = "Initial distribution")

s = 0
sw = FALSE
QQ = -Inf
QQ.out = NULL
epsilon = 10^(-5)

while (!sw) {
  v = array(0, dim = c(n, KK))
  v[, 1] = log(w) + log((x == 0) + 1e-3)
  v[, 2] = log(1 - w) + dpois(x, lambda, log = TRUE)
  for (i in 1:n) {
    v[i, ] = exp(v[i, ] - max(v[i, ])) / sum(exp(v[i, ] - max(v[i, ])))
  }

  w = mean(v[, 1])

  lambda = 0
  for (i in 1:n) {
    lambda = lambda + v[i, 2] * x[i]
  }
  lambda = lambda / sum(v[, 2])

  QQn = 0
  for (i in 1:n) {
    QQn = QQn + v[i, 1] * (log(w)) + v[i, 2] * (log(1 - w) + dpois(x[i], lambda, log = TRUE))
  }
  if (abs(QQn - QQ) / abs(QQn) < epsilon) {
    sw = TRUE
  }
  QQ = QQn
  QQ.out = c(QQ.out, QQ)
  s = s + 1
  print(paste(s, QQn))

  par(mfrow = c(3, 1))
  plot(QQ.out[1:s], type = "l", xlim = c(1, max(10, s)), las = 1, ylab = "Q", lwd = 2)

  xx = seq(0, max(x))
  yy = w * (xx == 0) + (1 - w) * dpois(xx, lambda)
  plot(xx, yy, type = "h", ylim = c(0, max(yy)), main = paste("s =", s, "   Q =", round(QQ.out[s], 4)), col = "red", xlab = "x", ylab = "Density")
  plot(table(x) / n, xlab = "x", ylab = "Density")
}

par(mfrow = c(3, 1))
plot(QQ.out[1:s], type = "l", xlim = c(1, max(10, s)), las = 1, ylab = "Q", lwd = 2)

xx = seq(0, max(x))
yy = w * (xx == 0) + (1 - w) * dpois(xx, lambda)
plot(xx, yy, type = "h", ylim = c(0, max(yy)), main = paste("s =", s, "   Q =", round(QQ.out[s], 4)), col = "red", xlab = "x", ylab = "Density")
plot(table(x) / n, xlab = "x", ylab = "Density")
