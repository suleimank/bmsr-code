rm(list=ls())

library(rstan)
source("runFuncs.R")
set.seed(101)


S = 2; nY = sample(20:30,S,replace = T);
dX = 50; dY = 3;
ft = 10; sigma = seq(0.1,0.3,l=S)
X = matrix(rnorm(sum(nY)*dX),sum(nY),dX)
betaM = rep(0,dX); betaM[1:ft] = (rnorm(ft,0,1))

#tau = sample(c(1,1e-3),S,replace = T)
Beta = matrix(0,S,dX); Y = matrix(NA,sum(nY),dY); W = matrix(runif(dY,0.8,1),1,dY)
for( s in 1:S){
  Beta[s,1:ft] = abs(rnorm(ft,0,1)); Beta[s,1:ft] = Beta[s,1:ft]/max(Beta[s,1:ft])
  Beta[s,1:ft] = Beta[s,1:ft]*betaM[1:ft]#*tau[s]
  if(s == 1) { st=1; en=nY[1]; }
  else { st=sum(nY[1:(s-1)])+1; en=sum(nY[1:s]); }
  Y[st:en,] =  (X[st:en,] %*% Beta[s,] %*% W);
  Y[st:en,] = Y[st:en,] + matrix(apply(Y[st:en,],2,sd),nY[s],dY,byrow=T)*sigma[s]*matrix(rnorm(nY[s]*dY,0,1),nY[s],dY)
}

data = list(S = S, nY = nY,dX = dX, dY = dY, Y = Y, X = X, p0 = ft*1.2)

opts = list(iter=300,seeds=c(12,345,6789),inference="Sampling")

file = "bmsr.stan";
predFunction = paste0("predict.",file);
posteriorFunction = paste0("posterior.",file);

xTest=X; nTest = nY 

res = runSTAN(file,data,opts)
out = res$out
yPred = predictSTAN(predFunction,out,xTest,nTest,yN)
print(diag(cor(yPred,Y)))
print( sqrt(mean( (yPred - Y)^2 )) )
post = getPosterior(posteriorFunction,out)

tt = c(as.vector(post$beta),as.vector(Beta))
op <- par(mfrow=c(S,1),mgp=c(0,0,0),mar=c(0,0.5,0,0)+0.5)
for(s in 1:S){
  plot(Beta[s,],ylim=c(min(tt),max(tt)),pch=16,xlab="",ylab="",tck=0.01)
  points(post$beta[s,],col="red",pch=17)
}
par(op)

