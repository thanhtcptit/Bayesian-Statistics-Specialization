n = 200
w = c(0.7, 0.2, 0.1)
mean = c(1, 2, 6)
cc = sample(1:3, n, replace = T, prob = w)

x = rpois(n, mean[cc])
hist(x)
