source("HelpingFunctions.R")
require("rstan")

#' train a regression STAN model
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
  #ml = try(load(file))
  #if(class(ml)== "try-error"){
  print("Compiling Stan Model")
  model = rstan::stan_model(file=list.files(pattern=file))
  #}
  
  print("Training Stan Model")
  ptm <- proc.time()
  if(opts$iter < 1e3)
    print(paste('Iterations >= 1e3 recommended, currently set to:',opts$iter))
  sampling_iterations = opts$iter
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

############################ PREDICT FUNCTIONS #########################################

#' predict from a regression STAN model
#'
#' \code{predictSTAN} predicts the output a stan regression model as defined by the parameter \code{predFunction}.
#'
#'
#' @param predFunction is a function for predicting y's using the outcome of stan run.
#' @param xTest is a matrix of test data for predicting the outcome. If NULL no prediction is made (default).
#' @param nTest is a vector of length S, containing the number of values in each source. Can contain zero's.
#' @param yN a list containing values used for normalizing the data: (default = NULL)
#'  \item{cm} vector of means with which the data is centered. 0's if data is not centered
#'  \item{cs} vector of standard deviations with which the data is scaled. 1's if data is not scaled.
#' @return yPred prediction vector of the stan model.
#' @export
predictSTAN <- function(predFunction,out,xTest,nTest,yN=NULL){
  yPred = NA;
  if(length(xTest)>0) {
    yPred = get(predFunction)(out, xTest, nTest, yN)
  }
  return(yPred)
}

#' predict function for single source
#'
#'
#' @param xTest.S is a matrix of test data for predicting the outcome.
#' @param beta is model parameter vector.
#' @param W are model parameters. If NULL, model is single/task, else multi/task.
#' @param yN a list containing values used for normalizing the data: (default = NULL)
#'  \item{cm} vector of means with which the data is centered. 0's if data is not centered
#'  \item{cs} vector of standard deviations with which the data is scaled. 1's if data is not scaled.
#' @return yPred prediction vector of the stan model.
pred <- function(xTest.S,beta,W,yN)
{
  if(length(W)>0)
    ynew = (xTest.S %*% t(beta) %*% W / nrow(W))
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

#' predict function for bmsr model
#'
#' \code{predict.bmsr.stan} predicts the output of bmsr model, used as \code{predFunction} in \code{predictSTAN}.
#'
#'
#' @param out is trained STAN model.
#' @param xTest is a matrix of test data for predicting the outcome. If NULL no prediction is made (default).
#' @param nTest is a vector of length S, containing the number of values in each source. Can contain zero's.
#' @param yN a list containing values used for normalizing the data: (default = NULL)
#'  \item{cm} vector of means with which the data is centered. 0's if data is not centered
#'  \item{cs} vector of standard deviations with which the data is scaled. 1's if data is not scaled.
#' @return yPred prediction vector of the bmsr model.
#' @export
predict.bmsr.stan <- function(out, xTest, nTest, yN=NULL)
{
  post = list()
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

#' predict function for bmsmtr model
#'
#' \code{predict.bmsmtr.stan} predicts the output of bmsmtr model, used as \code{predFunction} in \code{predictSTAN}.
#'
#'
#' @param out is trained STAN model.
#' @param xTest is a matrix of test data for predicting the outcome. If NULL no prediction is made (default).
#' @param nTest is a vector of length S, containing the number of values in each source. Can contain zero's.
#' @param yN a list containing values used for normalizing the data: (default = NULL)
#'  \item{cm} vector of means with which the data is centered. 0's if data is not centered
#'  \item{cs} vector of standard deviations with which the data is scaled. 1's if data is not scaled.
#' @return yPred prediction vector of the bmsmtr model. 
#' @export
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

############################################ GET POSTERIOR FUNCTIONS ##############################################

#' baseline function to get posterior
#'
#' \code{getPosterior} extracts the posterior values from mode output.
#'
#'
#' @param out is trained STAN model.
#' @param file is the stan file name containig the stan code.
#' @return post is a list containing posterior of all model weights.
getPosterior <- function(file=NULL,out)
{
  return(get(file)(out))
}

#'  get posterior of the bmsr model weights
#'
#' \code{posterior.bmsr.stan} extracts the posterior values from mode output.
#'
#' @param out is trained STAN model.
#' @return post is a list containing following posterior weights.
#'  \item{beta} matrix containing source specific beta parameters.
#'  \item{betaShared} vector of shared beta parameters common for all sources.
#'  \item{tau} global scaling factor learned.
#'  \item{sigma} noise parameter.
#' @export
posterior.bmsr.stan <- function(out)
{
  post = list()
  post$sigma = rstan::extract(out,"sigma",permuted=TRUE); 
  post$sigma = apply(post$sigma[[1]],c(2),mean)

  post$tau = rstan::extract(out,"tauM",permuted=TRUE);
  post$tau = mean(post$tau[[1]])
  
  post$beta = rstan::extract(out,"betaT",permuted=TRUE)
  post$beta = apply(post$beta[[1]],c(2,3),mean)

  tau0M = rstan::extract(out,"tau0M",permuted=TRUE)[[1]];
  tauM = rstan::extract(out,"tauM",permuted=TRUE)[[1]];  
  lambda = rstan::extract(out,"lambda",permuted=TRUE)[[1]];  
  betaM = rstan::extract(out,"betaM",permuted=TRUE)[[1]];
  post$betaShared <- t(lambda * betaM) %*% (tauM * tau0M)
  
  return(post)
}

#'  get posterior of the bmsmtr model weights
#'
#' \code{getPosterior} extracts the posterior values from mode output.
#'
#' @param out is trained STAN model.
#' @return post is a list containing following posterior weights.
#'  \item{beta} matrix containing source specific beta parameters.
#'  \item{W} matrix of multi-task parameters.
#'  \item{tau} global scaling factor learned.
#'  \item{sigma} noise parameter.
#' @export
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

#'  get beta posterior of the bmsr model
#'
#' \code{getBeta.bmsr.stan} extracts the posterior values of source specific beta parameters.
#'
#' @param out is trained STAN model.
#' @return beta matrix containing source specific beta parameters.
#' @export
getBeta.bmsr.stan <- function(out)
{
  beta = rstan::extract(out,"betaT",permuted=TRUE)
  beta = t(apply(beta[[1]],c(2,3),mean))
  return(beta)
}

#'  get beta posterior of the bmsmtr model
#'
#' \code{getBeta.bmsmtr.stan} extracts the posterior values of source specific beta parameters.
#'
#' @param out is trained STAN model.
#' @return beta matrix containing source specific beta parameters.
#' @export
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