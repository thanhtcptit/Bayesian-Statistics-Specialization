library(mvtnorm)
library(MCMCpack)

dat = load("dataset/banknoteclassification.Rdata")

n = dim(banknote.training)[1]
m = dim(banknote.test)[1]
x = rbind(as.matrix(banknote.training), as.matrix(banknote.test))
y_train = as.numeric(banknote.training.labels == "genuine") + 1
y_test  = as.numeric(banknote.test.labels == "genuine") + 1

KK = 2
p = dim(x)[2]
w = rep(1, KK) / KK
mu = rmvnorm(KK, apply(x, 2, mean), var(x))
Sigma = array(0, dim = c(KK, p, p))
Sigma[1, , ] = var(x) / KK
Sigma[2, , ] = var(x) / KK
Sigma[3, , ] = var(x) / KK
cc = c(y_train, sample(1:KK, m, replace = TRUE, prob = w))

aa = rep(1, KK)
dd = apply(x, 2, mean)
DD = 10 * var(x)
nu = p
SS = var(x) / 3

burn = 100
rrr = 1000

cc.out = array(0, dim = c(rrr, n + m))
w.out = array(0, dim = c(rrr, KK))
mu.out = array(0, dim = c(rrr, KK, p))
Sigma.out = array(0, dim = c(rrr, KK, p, p))
logpost = rep(0, rrr)

for (s in 1:rrr) {
  for (i in (n + 1):(n + m)) {
    v = rep(0, KK)
    for (k in 1:KK) {
      v[k] = log(w[k]) + mvtnorm::dmvnorm(x[i, ], mu[k, ], Sigma[k, , ], log = TRUE)
    }
    v = exp(v - max(v)) / sum(exp(v - max(v)))
    cc[i] = sample(1:KK, 1, replace = TRUE, prob = v)
  }

  w = as.vector(rdirichlet(1, aa + tabulate(cc)))

  DD.st = matrix(0, nrow=p, ncol=p)
  for (k in 1:KK) {
    mk = sum(cc == k)
    xsumk = apply(x[cc == k, ], 2, sum)
    DD.st = solve(mk * solve(Sigma[k, , ]) + solve(DD))
    dd.st = DD.st %*% (solve(Sigma[k, , ]) %*% xsumk + solve(DD) %*% dd)
    mu[k, ] = as.vector(rmvnorm(1, dd.st, DD.st))
  }

  xcensumk = array(0, dim = c(KK, p, p))
  for (i in 1:(n + m)) {
    xcensumk[cc[i], , ] = xcensumk[cc[i], , ] + (x[i, ] - mu[cc[i], ]) %*% t(x[i, ] - mu[cc[i], ])
  }
  for (k in 1:KK) {
    Sigma[k, , ] = riwish(nu + sum(cc == k), SS + xcensumk[k, , ])
  }

  cc.out[s, ] = cc
  w.out[s, ] = w
  mu.out[s, , ] = mu
  Sigma.out[s, , , ] = Sigma
  for (i in 1:(n + m)) {
    logpost[s] = logpost[s] + log(w[cc[i]]) + mvtnorm::dmvnorm(x[i, ], mu[cc[i], ], Sigma[cc[i], , ], log = TRUE)
  }
  logpost[s] = logpost[s] + ddirichlet(w, aa)
  for (k in 1:KK) {
    logpost[s] = logpost[s] + mvtnorm::dmvnorm(mu[k, ], dd, DD, log = TRUE)
    logpost[s] = logpost[s] + log(diwish(Sigma[k, , ], nu, SS))
  }

  if (s / 250 == floor(s / 250)) {
    print(paste("s = ", s))
  }
}

y_pred = as.integer(apply(cc.out[-(1:burn), ], 2, mean)[(n + 1):(n + m)])
sum(!(y_pred == y_test))

modqda = qda(grouping=y_train, x=banknote.training, method = "mle")
ccpredqda = predict(modqda, newdata = banknote.test)
sum(!(ccpredqda$class == y_test))