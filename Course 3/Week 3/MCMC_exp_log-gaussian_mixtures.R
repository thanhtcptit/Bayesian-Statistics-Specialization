rm(list = ls())
library(MCMCpack)
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
tau = sd(log(x))

xx = seq(0, 10, length = 200)
yy = w * dexp(xx, lambda) + (1 - w) * dlnorm(xx, mu, tau)
plot(xx, yy, type = "l", ylim = c(0, max(yy)), xlab = "x", ylab = "Initial density")
points(x, rep(0, n))

rrr = 6000
burn = 1000

lambda.p = 1
mu.p_mu = 0
mu.p_sigma2 = 1
tau.p_alpha = 2
tau.p_beta = 1

cc.out = array(0, dim = c(rrr, n))
w.out = rep(0, rrr)
lambda.out = rep(0, rrr)
mu.out = rep(0, rrr)
tau.out = rep(0, rrr)
logpost = rep(0, rrr)

for (s in 1:rrr) {
  cc = rep(0, n)
  for (i in 1:n) {
    v = rep(0, KK)
    v[1] = log(w) + dexp(x[i], lambda, log = TRUE)
    v[2] = log(1 - w) + dlnorm(x[i], mu, tau, log = TRUE)
    v = exp(v - max(v)) / sum(exp(v - max(v)))
    cc[i] = sample(1:KK, 1, replace = TRUE, prob = v)
  }

  n1 = sum(cc == 1)
  n2 = sum(cc == 2)
  w = rbeta(1, n1 + 1, n2 + 1)
  lambda = rgamma(1, n1 + 1, sum(x[cc == 1]) + lambda.p)

  sigma2_star = 1 / (n2 / tau ^ 2 + 1 / mu.p_sigma2)
  mu_star = sigma2_star * (sum(log(x[cc == 2])) / tau ^ 2 + mu.p_mu / mu.p_sigma2)
  mu = rnorm(1, mu_star, sqrt(sigma2_star))

  tau = sqrt(rinvgamma(1, n2 / 2 + tau.p_alpha, 1 / 2 * sum((log(x[cc == 2]) - mu) ^ 2) + tau.p_beta))

  cc.out[s, ] = cc
  w.out[s] = w
  lambda.out[s] = lambda
  mu.out[s] = mu
  tau.out[s] = tau
  for (i in 1:n) {
    if (cc[i] == 1) {
      logpost[s] = logpost[s] + log(w) + dexp(x[i], lambda, log = TRUE)
    } else {
      logpost[s] = logpost[s] + log(1 - w) + dlnorm(x[i], mu, tau, log = TRUE)
    }
  }
  logpost[s] = logpost[s] + dexp(lambda, lambda.p, log = TRUE)
  logpost[s] = logpost[s] + dnorm(mu, mu.p_mu, sqrt(mu.p_sigma2), log = TRUE)
  logpost[s] = logpost[s] + log(dinvgamma(tau^2, tau.p_alpha, tau.p_beta))
  if (s / 500 == floor(s / 500)) {
    print(paste("s =", s))
  }
}

mean_w = mean(w.out[-seq(1:burn)])
mean_lambda = mean(lambda.out[-seq(1:burn)])
mean_mu = mean(mu.out[-seq(1:burn)])
mean_tau = mean(tau.out[-seq(1:burn)])

xx = seq(0, 10, length = 200)
yy = mean_w * dexp(xx, mean_lambda) + (1 - mean_w) * dlnorm(xx, mean_mu, mean_tau)
plot(xx, yy, type = "l", ylim = c(0, max(yy)), col = "red", xlab = "x", ylab = "Density")
points(x, rep(0, n))
