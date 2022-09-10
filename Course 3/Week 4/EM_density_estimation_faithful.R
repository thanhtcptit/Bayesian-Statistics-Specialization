rm(list = ls())

library(MASS)
library(MCMCpack)
data(faithful)

x = faithful$eruptions
n = length(x)
set.seed(781209)

KK = 2
w = rep(1, KK) / KK
mu = rnorm(KK, mean(x), sd(x))
sigma = c(sd(x) / KK, sd(x) / KK)

epsilon = 0.000001
s = 0
sw = FALSE
KL = -Inf
KL.out = NULL

while (!sw) {
  v = array(0, dim = c(n, KK))
  for (k in 1:KK) {
    v[, k] = log(w[k]) + dnorm(x, mu[k], sigma[k], log = TRUE)
  }
  for (i in 1:n) {
    v[i, ] = exp(v[i, ] - max(v[i, ])) / sum(exp(v[i, ] - max(v[i, ])))
  }

  w = apply(v, 2, mean)
  mu = rep(0, KK)
  for (k in 1:KK) {
    for (i in 1:n) {
      mu[k] = mu[k] + v[i, k] * x[i]
    }
    mu[k] = mu[k] / sum(v[, k])
  }

  for (k in 1:KK) {
    sigma[k] = 0
    for (i in 1:n) {
        sigma[k] = sigma[k] + v[i, k] * (x[i] - mu[k]) ^ 2
    }
    sigma[k] = sqrt(sigma[k] / sum(v[, k]))
  }

  ## Check convergence
  KLn = 0
  for (i in 1:n) {
    for (k in 1:KK) {
      KLn = KLn + v[i, k] * (log(w[k]) + dnorm(x[i], mu[k], sigma[k], log = TRUE))
    }
  }
  if (abs(KLn - KL) / abs(KLn) < epsilon) {
    sw = TRUE
  }
  KL = KLn
  KL.out = c(KL.out, KL)
  s = s + 1
  print(paste(s, KLn))
}
xx = seq(0, 7, length = 150)
nxx = length(xx)
density.EM = rep(0, nxx)
for (s in 1:nxx) {
  for (k in 1:KK) {
    density.EM[s] = density.EM[s] + w[k] * dnorm(xx[s], mu[k], sigma[k])
  }
}
plot(xx, density.EM, type="l")
points(x, rep(0, n))
