rm(list = ls())
set.seed(81196)

dat = read.csv("../dataset/fuses.csv", header = FALSE)
x = dat$V1
x_sort = sort(x)
n = length(x)
plot(x_sort, type = "h")

KK = 2
w = 1 / 2
lambda = 1 / mean(x)
mu = mean(log(x))
sigma = sd(log(x))

xx = seq(0, 10, length = 200)
yy = w * dexp(xx, lambda) + (1 - w) * dlnorm(xx, mu, sigma)
plot(xx, yy, type = "l", ylim = c(0, max(yy)), xlab = "x", ylab = "Initial density")
points(x, rep(0, n))

s = 0
sw = FALSE
QQ = -Inf
QQ.out = NULL
epsilon = 10^(-5)

while (!sw) {
  v = array(0, dim = c(n, KK))
  v[, 1] = log(w) + dexp(x, lambda, log = TRUE)
  v[, 2] = log(1 - w) + dlnorm(x, mu, sigma, log = TRUE)
  for (i in 1:n) {
    v[i, ] = exp(v[i, ] - max(v[i, ])) / sum(exp(v[i, ] - max(v[i, ])))
  }

  w = mean(v[, 1])

  lambda = 0
  for (i in 1:n) {
    lambda = lambda + v[i, 1] * x[i]
  }
  lambda = sum(v[, 1]) / lambda

  mu = 0
  for (i in 1:n) {
    mu = mu + v[i, 2] * log(x[i])
  }
  mu = mu / sum(v[, 2])

  sigma = 0
  for (i in 1:n) {
    sigma = sigma + v[i, 2] * (log(x[i]) - mu)^2
  }
  sigma = sqrt(sigma / sum(v[, 2]))

  QQn = 0
  for (i in 1:n) {
    QQn = QQn + v[i, 1] * (log(w) + dexp(x[i], lambda, log = TRUE)) +
      v[i, 2] * (log(1 - w) + dlnorm(x[i], mu, sigma, log = TRUE))
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

  xx = seq(0, 10, length = 200)
  yy = w * dexp(xx, lambda) + (1 - w) * dlnorm(xx, mu, sigma)
  plot(xx, yy, type = "l", ylim = c(0, max(yy)), main = paste("s =", s, "   Q =", round(QQ.out[s], 4)), col = "red", xlab = "x", ylab = "Density")
  points(x, rep(0, n))
}

par(mfrow = c(3, 1))
plot(QQ.out[1:s], type = "l", xlim = c(1, max(10, s)), las = 1, ylab = "Q", lwd = 2)

xx = seq(0, 10, length = 200)
yy = w * dexp(xx, lambda) + (1 - w) * dlnorm(xx, mu, sigma)
plot(xx, yy, type = "l", ylim = c(0, max(yy)), main = paste("s =", s, "   Q =", round(QQ.out[s], 4)), col = "red", xlab = "x", ylab = "Density")
points(x, rep(0, n))
