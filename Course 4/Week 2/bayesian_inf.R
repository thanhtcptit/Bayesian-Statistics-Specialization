# Simulate 300 observations from an AR(2) with one pair of complex-valued roots 
set.seed(2021)
r=0.95
lambda=12 
phi=numeric(2) 
phi[1]=2*r*cos(2*pi/lambda) 
phi[2]=-r^2
sd=1 # innovation standard deviation
T=300 # number of time points
# generate stationary AR(2) process
yt=arima.sim(n = T, model = list(ar = phi), sd = sd) 
par(mfrow=c(1,1))
plot(yt)

## Compute the MLE of phi and the unbiased estimator of v using the conditional likelihood
p=2
y=rev(yt[(p+1):T]) # response
X=t(matrix(yt[rev(rep((1:p),T-p)+rep((0:(T-p-1)),rep(p,T-p)))],p,T-p));
XtX=t(X)%*%X
XtX_inv=solve(XtX)
phi_MLE=XtX_inv%*%t(X)%*%y # MLE for phi
s2=sum((y - X%*%phi_MLE)^2)/(length(y) - p) #unbiased estimate for v

#####################################################################################
### Posterior inference, conditional likelihood + reference prior via 
### direct sampling                 
#####################################################################################

n_sample=1000 # posterior sample size
library(MASS)

## step 1: sample v from inverse gamma distribution
v_sample=1/rgamma(n_sample, (T-2*p)/2, sum((y-X%*%phi_MLE)^2)/2)

## step 2: sample phi conditional on v from normal distribution
phi_sample=matrix(0, nrow = n_sample, ncol = p)
for(i in 1:n_sample){
  phi_sample[i, ]=mvrnorm(1,phi_MLE,Sigma=v_sample[i]*XtX_inv)
}

par(mfrow = c(2, 3), cex.lab = 1.3)
## plot histogram of posterior samples of phi and v

for(i in 1:2){
  hist(phi_sample[, i], xlab = bquote(phi), 
       main = bquote("Histogram of "~phi[.(i)]),col='lightblue')
  abline(v = phi[i], col = 'red')
}

hist(v_sample, xlab = bquote(nu), main = bquote("Histogram of "~v),col='lightblue')
abline(v = sd, col = 'red')

#####################################################
# Graph posterior for modulus and period 
#####################################################
r_sample=sqrt(-phi_sample[,2])
lambda_sample=2*pi/acos(phi_sample[,1]/(2*r_sample))
hist(r_sample,xlab="modulus",main="",col='lightblue')
abline(v=0.95,col='red')
hist(lambda_sample,xlab="period",main="",col='lightblue')
abline(v=12,col='red')


