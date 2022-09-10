library("rjags")

# Load & Preprocess data
dat = read.csv("cpu.csv", header=FALSE)
colnames(dat) = c("vendor", "model", "myct", "mmin", "mmax", "cache", "chmin", "chmax", "prp", "erp")

colnames_subset = c("vendor", "myct", "mmax", "chmax", "cache", "prp")
d = dat[colnames_subset]
summary(d)

d$vendor = unclass(as.factor(d$vendor))

pairs(d)
par(mfrow=c(2, 3))
hist(d$myct, main = "myct")
hist(d$mmax, main = "mmax")
hist(d$chmax, main = "chmax")
hist(d$cache, main = "cache")
hist(d$prp, main = "prp")

log_transform_cols = c("myct", "mmax", "chmax", "cache", "prp")
for (i in log_transform_cols) {
    d[[paste("log", i, sep="")]] = log(d[[i]] + 1.00001)
}
d$ilogmyct = 1 / log(d$myct)
d$ilogmyct_logmmax_logchmax = d$ilogmyct * d$logmmax * d$logchmax

pairs(d[c("ilogmyct_logmmax_logchmax", "logcache", "logprp")])
par(mfrow=c(1, 3))
hist(d$ilogmyct_logmmax_logchmax, main = "ilogmyct_logmmax_logchmax")
hist(d$logcache, main = "logcache")
hist(d$logprp, main = "logprp")


# Define JAGS model
set.seed(42)
data1_jags = as.list(d)

mod1_string = " model {
    for (i in 1:length(logprp)) {
        logprp[i] ~ dnorm(mu[i], prec)
        mu[i] = a[vendor[i]] + b[1] * ilogmyct_logmmax_logchmax[i] + b[2] * logcache[i]
    }

    for (i in 1:max(vendor)) {
        a[i] ~ dnorm(0.0, 1.0 / 1.0e6)
    }
    for (i in 1:2) {
        b[i] ~ dnorm(0.5, 1.0 / 1.0e6)
    }
    prec ~ dgamma(5 / 2.0, 5 * 10.0 / 2.0)
    sig = sqrt(1 / prec)

} "

params1 = c("a", "b", "sig")
inits1 = function() {
    inits = list("a"=rnorm(max(d$vendor), 0.0, 100.0), "b"=rnorm(2, 0.5, 100.0), "prec"=rgamma(1, 1.0, 1.0))
}

mod1 = jags.model(textConnection(mod1_string), data=data1_jags, inits=inits1, n.chains=3)
update(mod1, 5000)

mod1_sim = coda.samples(model=mod1, variable.names=params1, n.iter=1e4)
mod1_csim = do.call(rbind, mod1_sim)


# Check coverage
plot(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
gelman.diag(mod1_sim)
summary(mod1_sim)


# Visualize residual
pm_params1 = colMeans(mod1_csim)
X = cbind(data1_jags$ilogmyct_logmmax_logchmax, data1_jags$logcache)
yhat1 = drop(X %*% pm_params1[31:32]) + pm_params1[data1_jags$vendor]
resid1 = data1_jags$logprp - yhat1

par(mfrow=c(1, 2))
plot(resid1)
plot(yhat1, resid1)
