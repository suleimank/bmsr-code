// Bayesian Multi-source sparse regression
//
// Author: Suleiman Ali Khan (c) 2017
//
// khan.suleiman@gmail.com
//
// Last Update: 28th Mar 2018

data {
  int<lower=0> S;
  int<lower=0> nY[S];
  int<lower=0> dX;
  int<lower=0> dY;
  matrix[sum(nY),dY] Y;
  matrix[sum(nY),dX] X;
  int<lower=0> p0;
}
transformed data {
  int<lower=0> N;
  int<lower=0> iY[S,2];
    
  N = sum(nY);
  for(s in 1:S){
    if(s == 1){
      iY[s,1] = 1;
      iY[s,2] = nY[s];
    }
    else{
      iY[s,1] = iY[s-1,2]+1;
      iY[s,2] = iY[s-1,2]+nY[s];
    }
  }
}

parameters {
  vector[dX] betaM;
  real<lower=0> tauM;
  
  matrix[S,dX] beta;
  vector<lower=0>[dX] lambda; 
  real<lower=0> sigma[S];
}

transformed parameters {
  vector[dY] W;
  real<lower=0> tau0M;
  matrix[S,dX] betaT;  
  
  for (d in 1:dY)  W[d] = 1;
  
  tau0M = (sum(sigma)*p0)/((dX-p0)*sqrt(N*1.0));

  for (s in 1:S) {
     betaT[s] = beta[s] .* (betaM .* lambda * tauM * tau0M)' ;
  }
}

model{
  lambda ~ cauchy(0,1);
  sigma ~ normal(0,1); // 18.12.2017 Previously used: cauchy(0,2.5);
  
  tauM ~ cauchy(0,1); 
  betaM ~ normal(0,1); 
  
  for (s in 1:S){
    beta[s] ~ normal(0.5,0.5);
    for(d in 1:dY){
      Y[iY[s,1]:iY[s,2],d] ~ normal(X[iY[s,1]:iY[s,2]]*(betaT[s]')*W[d],sigma[s]);
    }
  }
}
