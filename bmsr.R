source("HelpingFunctions.R")
require("rstan")

#' run a regression STAN model
#'
#' \code{runSTAN} runs a stan regression model and predicts the values for the test samples.
#'
#'
#' @param file is the stan file containig the stan code.
#' @param data is a list containing the data in the format this stan code accepts
#' @param opts a list containing opts to run the model:
#'  \item{iter} integer identifies number of sampling iterations
#'  \item{seeds} vector of integers identifying seeds for running the method
#'  \item{inference} string identifing the sampling method to use, either of Sampling or VB
#' @return A list containing the following elements:
#'  \item{out} STAN output variable
#'  \item{runtime} run time of the code. 
#' @export
runSTAN <- function(file,data,opts)
{
  ml = try(load(file))
  if(class(ml)== "try-error"){
    print("Compiling Stan...")
    model = rstan::stan_model(file=list.files(pattern=file))
  }
  
  print("Running Stan...")
  ptm <- proc.time()
  sampling_iterations = opts$iter #best to use 1e3 or higher
  seeds = opts$seeds 
  if(opts$inference == "Sampling")
  {
    for(i in 1:length(seeds))
    {
      out = try(rstan::sampling(
        object = model, data = data, chains = 1
        , iter = sampling_iterations, warmup = sampling_iterations/2, refresh = sampling_iterations/100 #show an update @ each %10
        , control = list(adapt_delta = 0.999, max_treedepth = 20), seed = seeds[i]))
      if(class(out)!="try-error") break;
    }
  }
  if(opts$inference == "VB")
  {
    for(i in 1:length(seeds))
    {
      out = try(vb(model, data = data, algorithm = "meanfield",tol_rel_obj = 1e-3,iter = sampling_iterations,seed=seeds[i],eval_elbo=1000))
      if(class(out)!="try-error") break;
    }
  }
  
  rt = proc.time() - ptm
  print(rt)
  return(list(out=out,runtime=rt))
}

###################################################################################################################################

#' @param predFunction is a function for predicting y's using the outcome of stan run.
#' @param xTest is a matrix of test data for predicting the outcome. If NULL no prediction is made (default).
#' @param nTest is a vector of length S, containing the number of values in each source. Can contain zero's.
#' @param yN a list containing values used for normalizing the data: (default = NULL)
#'  \item{cm} vector of means with which the data is centered. 0's if data is not centered
#'  \item{cs} vector of standard deviations with which the data is scaled. 1's if data is not scaled.
#' @export
predictSTAN <- function(predFunction,out,xTest,nTest,yN=NULL){
  yPred = NA;
  if(length(xTest)>0) {
    yPred = get(predFunction)(out, xTest, nTest, yN)
  }
  return(yPred)
}

pred <- function(xTest.S,beta,W,yN)
{
  if(length(W)>0)
    ynew = (xTest.S %*% t(beta) %*% W / nrow(W)) #same as for loop below
  else
    ynew = apply(xTest.S %*% t(beta),1,mean)
  
  if(length(yN)>0){
    if(nrow(xTest.S)==1) {
      ynew = ynew*yN$cs + yN$cm;
    } else {
      ynew = ynew*matrix(yN$cs,nrow(ynew),ncol(ynew),byrow=T) + matrix(yN$cm,nrow(ynew),ncol(ynew),byrow=T);
    }
  }
  return(ynew)
}

predict.bmsr.stan <- function(out, xTest, nTest, yN=NULL)
{
  post = list()
  #post$W = rstan::extract(out,"W",permuted=TRUE)
  post$beta = rstan::extract(out,"betaT",permuted=TRUE)
  S = length(nTest)
  yNew = matrix(NA,nrow(xTest),1)
  for(s in 1:S){
    if(nTest[s]>0){
      if(s == 1) { st = 1; en = nTest[s]; }
      else { st = sum(nTest[1:(s-1)])+1 ; en = sum(nTest[1:s]); }
      yNew[st:en,] = pred(xTest[st:en,],post$beta[[1]][,s,],NULL,yN=NULL)
    }
  }
  return(yNew)
}

predict.bmsmtr.stan <- function(out, xTest, nTest, yN=NULL)
{
  post = list()
  post$W = rstan::extract(out,"W",permuted=TRUE)
  post$W[[1]] = (post$W[[1]]/2 + 0.5)
  print(dim(post$W[[1]]))
  post$beta = rstan::extract(out,"betaT",permuted=TRUE)
  S = length(nTest)
  yNew = matrix(NA,nrow(xTest),dim(post$W[[1]])[2])
  for(s in 1:S){
    if(nTest[s]>0){
      if(s == 1) { st = 1; en = nTest[s]; }
      else { st = sum(nTest[1:(s-1)])+1 ; en = sum(nTest[1:s]); }
      yNew[st:en,] = pred(xTest[st:en,],post$beta[[1]][,s,],post$W[[1]],yN=NULL)
    }
  }
  return(yNew)
}

#######################################################################################################################################

getPosterior <- function(file=NULL,out)
{
  return(get(file)(out))
}

posterior.bmsr.stan <- function(out)
{
  post = list()
  post$sigma = rstan::extract(out,"sigma",permuted=TRUE); 
  post$sigma = apply(post$sigma[[1]],c(2),mean)
  
  post$tau = rstan::extract(out,"tauM",permuted=TRUE);
  post$tau = mean(post$tau[[1]])
  
  #post$W = rstan::extract(out,"W",permuted=TRUE);
  #post$W = apply(post$W[[1]],c(2),mean)
  #post$W = matrix(post$W,1,length(post$W))
  
  post$beta = rstan::extract(out,"betaT",permuted=TRUE)
  post$beta = apply(post$beta[[1]],c(2,3),mean)
  
  #sc = sum(post$W);
  #post$W = post$W * 1./sc
  #post$beta = post$beta * sc
  
  return(post)
}

posterior.bmsmtr.stan <- function(out)
{
  post = list()
  post$sigma = rstan::extract(out,"sigma",permuted=TRUE); 
  post$sigma = apply(post$sigma[[1]],c(2),mean)
  
  post$tau = rstan::extract(out,"tauM",permuted=TRUE);
  post$tau = mean(post$tau[[1]])
  
  post$W = rstan::extract(out,"W",permuted=TRUE);
  post$W = apply(post$W[[1]],c(2),mean)
  post$W = matrix(post$W,1,length(post$W))
  
  post$beta = rstan::extract(out,"betaT",permuted=TRUE)
  post$beta = apply(post$beta[[1]],c(2,3),mean)
  
  sc = sum(post$W);
  post$W = post$W * 1./sc
  post$beta = post$beta * sc
  
  return(post)
}

getBeta.bmsr.stan <- function(out)
{
  beta = rstan::extract(out,"betaT",permuted=TRUE)
  beta = t(apply(beta[[1]],c(2,3),mean))
  return(beta)
}

getBeta.bmsmtr.stan <- function(out)
{
  beta = rstan::extract(out,"betaT",permuted=TRUE)
  beta = t(apply(beta[[1]],c(2,3),mean))
  
  W = rstan::extract(out,"W",permuted=TRUE);
  W = apply(W[[1]],c(2),mean)
  W = matrix(W,1,length(W))
  
  sc = sum(W);
  W = W * 1./sc
  beta = beta * sc
  
  return(beta)
}