rm(list = ls())
set.seed(81196)

dat = read.csv("nestsize.csv", header = FALSE)
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

lambda_p = 1
rrr = 6000
burn = 1000

cc = rep(0, n)
cc.out = array(0, dim = c(rrr, n))
w.out = rep(0, rrr)
lambda.out = rep(0, rrr)
logpost = rep(0, rrr)

for (s in 1:rrr) {
  cc = rep(0, n)
  for (i in 1:n) {
    v = rep(0, KK)
    v[1] = log(w) + log((x[i] == 0) + 1e-5)
    v[2] = log(1 - w) + dpois(x[i], lambda, log = TRUE)
    v = exp(v - max(v)) / sum(exp(v - max(v)))
    cc[i] = sample(1:KK, 1, replace = TRUE, prob = v)
  }

  n2 = sum(cc == 2)
  w = rbeta(1, sum(cc == 1) + 1, n2 + 1)
  lambda = rgamma(1, sum(x[cc == 2]) + 1, n2 + lambda_p)

  cc.out[s, ] = cc
  w.out[s] = w
  lambda.out[s] = lambda

  for (i in 1:n) {
    if (cc[i] == 1) {
      logpost[s] = logpost[s] + log(w) + log((x[i] == 0) + 1e-5)
    } else {
      logpost[s] = logpost[s] + log(1 - w) + dpois(x[i], lambda, log = TRUE)
    }
    logpost[s] = logpost[s] - x[i]
  }
  if (s / 500 == floor(s / 500)) {
    print(paste("s =", s))
  }
}

mean_w = mean(w.out[-seq(1:burn)])
mean_lambda = mean(lambda.out[-seq(1:burn)])

xx = seq(0, max(x))
yy = mean_w * (xx == 0) + (1 - mean_w) * dpois(xx, mean_lambda)
par(mfrow = c(2, 1))
plot(table(x), type = "h")
plot(xx, yy, type = "h", ylim = c(0, max(yy)), xlab = "x", ylab = "Distribution from mean posterior")
