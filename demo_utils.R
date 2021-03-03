#' generate multi-source toy data for regression
#'
#' \code{generateSyntheticData} simualtes random multi-source multi-task regression dataset.
#'
#'
#' @param S is a scaler representing the desired number of sources (default = 2).
#' @param nY is a vector representing the desired number of samples in each source.
#' @param dX is a scaler representing the desired number of dimensions in X (inputs).
#' @param dY is a scaler representing the desired number of dimensions in Y (outputs). dY > 1 refers to multi-task datasets.
#' @param ft is a scaler representing the desired number of active features in X.
#' @return A list containing the following elements:
#'  \item{data} a list containing the Y and X data matrices
#'  \item{Beta} the beta parameters used to generate the data 
#' @export
generateSyntheticData <- function(  S = 2, 
                              nY = c(40,80),
                              dX = 50,
                              dY = 1,
                              ft = 10 )
{
  sprintf('    Generating synthetic data from S:%i, dX:%i, dY:%i',S,dX,dY)
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

#' demo multi-source model training and interpretative plots
#'
#' \code{demo_bmsr} trains bmsr and bmsmtr on random multi-source regression dataset.
#'
#'
#' @param file is the stan file containig the stan code.
#' @param dY is a scaler representing the desired number of dimensions in Y (outputs). dY > 1 refers to multi-task datasets.
#' @export
demo_bmsr <- function(file = "bmsr.stan",dY = 1)
{
  #data
  ll = generateSyntheticData(S = 2, 
                           nY = c(40,80),
                           dX = 50,
                           dY = dY,
                           ft = 10); 
  data = ll$data; Beta = ll$Beta
  xTest = data$X 
  nTest = data$nY 
  Y = data$Y
  S = data$S
  dX = ncol(xTest)
  dY = ncol(Y)
  
  #parameters
  opts = list(iter=1000,seeds=c(12,345,6789),inference="Sampling")
  predFunction = paste0("predict.",file);
  posteriorFunction = paste0("posterior.",file);
  
  #training
  sprintf('    Training Model %s on S:%i, dX:%i, dY:%i',file,S,dX,dY)
  res = runSTAN(file,data,opts)
  out = res$out
  
  #prediction
  sprintf('    Predict from Model %s on S:%i, dX:%i, dY:%i',file,S,dX,dY)
  yPred = predictSTAN(predFunction,out,xTest,nTest,yN)
  sprintf('    Performance of Model %s on S:%i, dX:%i, dY:%i',file,S,dX,dY)
  sprintf('    Correlation: %.3f',diag(cor(yPred,Y)))
  sprintf('    RMSE: %.3f',sqrt(mean( (yPred - Y)^2 )) )
  
  #interpretive plots
  post = getPosterior(posteriorFunction,out)
  tt = c(as.vector(post$beta),as.vector(Beta))
  op <- par(mfrow=c(S,1),mgp=c(0.5,0,0),mar=c(1,2,0.5,0)+0.5)
  cols = c('black','red')
  pchs = c(16,17)
  for(s in 1:S){
    plot(Beta[s,],ylim=c(min(tt),max(tt)),pch=pchs[1],xlab="feature in dX",ylab="weight",tck=0.01,main = paste0('Feature weights in source ',s),col=cols[1],cex.lab=.5, cex.axis=.5, cex.main=.5, cex.sub=.5)
    lines(c(-100,100), c(0,0) ,col="gray"); lines(c(0,0), c(-100,100) ,col="gray")
    points(post$beta[s,],col=cols[2],pch=pchs[2])
    legend("topright", 
           legend = c("true", "predicted"), 
           col = cols, 
           pch = pchs, 
           bty = "n",
           cex = 0.7, 
           text.col = "black", 
           horiz = F)
  }
  par(op)
  sprintf('    COMPLETED Model %s on S:%i, dX:%i, dY:%i',file,S,dX,dY)
}
