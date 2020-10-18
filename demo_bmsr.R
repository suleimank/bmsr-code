rm(list=ls())

library(rstan)
source("bmsr.R")
set.seed(101)

generateToyData <- function()
{
  S = 2; 
  nY = c(40,80);
  dX = 50; 
  dY = 1;
  ft = 10; 
  
  sigma = seq(0.2,0.5,l=S)
  X = matrix(rnorm(sum(nY)*dX),sum(nY),dX)
  betaM = rep(0,dX); betaM[1:ft] = rnorm(ft,0.7,1)
  
  Beta = matrix(1e-3,S,dX); 
  Y = matrix(NA,sum(nY),dY);
  for( s in 1:S){
    Beta[s,1:ft] = abs(rnorm(ft,0.5,0.5)); Beta[s,1:ft] = Beta[s,1:ft]/max(Beta[s,1:ft])
    Beta[s,1:ft] = Beta[s,1:ft]*betaM[1:ft]
    if(s == 1) { st=1; en=nY[1]; }
    else { st=sum(nY[1:(s-1)])+1; en=sum(nY[1:s]); }
    Y[st:en,] =  (X[st:en,] %*% Beta[s,]);
    Y[st:en,] = Y[st:en,] + matrix(sd(Y[st:en,]),nY[s],dY,byrow=T)*sigma[s]*matrix(rnorm(nY[s]*dY,0,1),nY[s],dY)
  }
  
  data = list(S = S, nY = nY,dX = dX, dY = 1, Y = matrix(Y[,1],nrow(Y),1), X = X, p0 = ft*2)
  
  return(list(data=data,Beta=Beta))
}

ll = generateToyData(); data = ll$data; Beta = ll$Beta
xTest = data$X 
nTest = data$nY 
Y = data$Y
dim(Y)
S = data$S

opts = list(iter=1000,seeds=c(12,345,6789),inference="Sampling")

file = "bmsr.stan";
predFunction = paste0("predict.",file);
posteriorFunction = paste0("posterior.",file);

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

