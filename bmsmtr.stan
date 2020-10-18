// Bayesian multi-source multi-task sparse regression
// 
// Citation: To Appear
//
// Author: Suleiman Ali Khan (c) 2019
//
// Email: khan.suleiman@gmail.com

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
  vector[dY] W;
  vector[dX] betaM;
  real<lower=0> tauM;
  matrix[S,dX] beta;
  vector<lower=0>[dX] lambda; 
  real<lower=0> sigma[S];
}

transformed parameters {
  real<lower=0> tau0M;
  matrix[S,dX] betaT;  
  tau0M = (sum(sigma)*p0)/((dX-p0)*sqrt(N*1.0));
  for (s in 1:S) {
     // betaT[s] = beta[s] .* (betaM .* lambda * tauM * tau0M)' ;
     betaT[s] = (beta[s] ./2 + 0.5) .* (betaM .* lambda * tauM * tau0M)' ;
  }
}

model{
  W ~ normal(0,1);
  lambda ~ cauchy(0,1);
  sigma ~  inv_gamma(1,1);
  tauM ~ cauchy(0,1); 
  betaM ~ normal(0,1); 
  
  for (s in 1:S){
    // beta[s] ~ normal(0.5,0.5);
    beta[s] ~ normal(0,1);
    for(d in 1:dY){
      Y[iY[s,1]:iY[s,2],d] ~ normal(X[iY[s,1]:iY[s,2]]*(betaT[s]')*(0.5 + W[d]/2),sigma[s]);
    }
  }
}
