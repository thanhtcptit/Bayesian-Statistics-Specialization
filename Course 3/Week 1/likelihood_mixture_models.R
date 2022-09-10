n = 100
w = c(0.3, 0.25, 0.25, 0.2)
mean = c(1, 4, 7, 10)
cc = sample(1:4, n, replace = T, prob = w)

x = rexp(n, 1 / mean[cc])
mean(x)
var(x)
