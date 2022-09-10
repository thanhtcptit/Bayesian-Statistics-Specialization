rm(list = ls())
set.seed(81196)

dat = read.csv("../dataset/nestsize.csv", header = FALSE)
x = dat$V1
n = length(x)

KK = 2
w = 1 / 2
lambda = mean(x)

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
}

BIC_mixture = 0
for(i in 1: n) {
    BIC_mixture = BIC_mixture - 2 * log(w * as.numeric(x[i] == 0) + (1 - w) * dpois(x[i], lambda))
}
BIC_mixture = BIC_mixture + ((KK - 1) + 1) * log(n)

lambda_sp = mean(x)
BIC_sp = 0
for(i in 1: n) {
    BIC_sp = BIC_sp - 2 * log(dpois(x[i], lambda_sp))
}
BIC_sp = BIC_sp + log(n)
