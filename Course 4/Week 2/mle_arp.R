set.seed(2021)
# Simulate 300 observations from an AR(2) with one pair of complex-valued reciprocal roots 
r=0.95
lambda=12 
phi=numeric(2) 
phi[1]=2*r*cos(2*pi/lambda) 
phi[2]=-r^2
sd=1 # innovation standard deviation
T=300 # number of time points
# generate stationary AR(2) process
yt=arima.sim(n = T, model = list(ar = phi), sd = sd) 

## Compute the MLE for phi and the unbiased estimator for v using the conditional likelihood
p=2
y=rev(yt[(p+1):T]) # response
X=t(matrix(yt[rev(rep((1:p),T-p)+rep((0:(T-p-1)),rep(p,T-p)))],p,T-p));
XtX=t(X)%*%X
XtX_inv=solve(XtX)
phi_MLE=XtX_inv%*%t(X)%*%y # MLE for phi
s2=sum((y - X%*%phi_MLE)^2)/(length(y) - p) #unbiased estimate for v

cat("\n MLE of conditional likelihood for phi: ", phi_MLE, "\n",
    "Estimate for v: ", s2, "\n")

