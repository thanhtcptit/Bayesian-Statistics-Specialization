###################################################
# Simulate data from an AR(2)
###################################################
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

#############################################################################
######   compute AIC and BIC for different AR(p)s based on simulated data ###
#############################################################################
pmax=10 # the maximum of model order
Xall=t(matrix(yt[rev(rep((1:pmax),T-pmax)+rep((0:(T-pmax-1)),
              rep(pmax,T-pmax)))], pmax, T-pmax));
y=rev(yt[(pmax+1):T])
n_cond=length(y) # (number of total time points - the maximum of model order)

## compute MLE
my_MLE <- function(y, Xall, p){
  n=length(y)
  x=Xall[,1:p]
  a=solve(t(x) %*%x)
  a=(a + t(a))/2 # for numerical stability 
  b=a%*%t(x)%*%y # mle for ar coefficients
  r=y - x%*%b # residuals 
  nu=n - p # degrees freedom
  R=sum(r*r) # SSE
  s=R/nu #MSE
  return(list(b = b, s = s, R = R, nu = nu))
}


## function for AIC and BIC computation 
AIC_BIC <- function(y, Xall, p){
  ## number of time points
  n <- length(y)
  
  ## compute MLE
  tmp=my_MLE(y, Xall, p)
  
  ## retrieve results
  R=tmp$R
  
  ## compute likelihood
  likl= n*log(R)
  
  ## compute AIC and BIC
  aic =likl + 2*(p)
  bic =likl + log(n)*(p)
  return(list(aic = aic, bic = bic))
}
# Compute AIC, BIC 
aic =numeric(pmax)
bic =numeric(pmax)

for(p in 1:pmax){
  tmp =AIC_BIC(y,Xall, p)
  aic[p] =tmp$aic
  bic[p] =tmp$bic
  print(c(p, aic[p], bic[p])) # print AIC and BIC by model order
}

## compute difference between the value and its minimum
aic =aic-min(aic) 
bic =bic-min(bic) 

## draw plot of AIC, BIC, and the marginal likelihood
par(mfrow = c(1, 1))
matplot(1:pmax,matrix(c(aic,bic),pmax,2),ylab='value',
        xlab='AR order p',pch="ab", col = 'black', main = "AIC and BIC")
# highlight the model order selected by AIC
text(which.min(aic), aic[which.min(aic)], "a", col = 'red') 
# highlight the model order selected by BIC
text(which.min(bic), bic[which.min(bic)], "b", col = 'red') 

########################################################
p <- which.min(bic) # We set up the moder order
print(paste0("The chosen model order by BIC: ", p))
