
// defining data
data {
  int<lower=1> N;
  int<lower=1> K; // number of groups
  real x[N]; 
  real y[N];
  int<lower=1,upper=K> Groups[N];
}

parameters {
  real<lower=0> slope[K];
  real<lower=0> sigma[K]; 
}

model{
for (i in 1:N)
  {
    target += normal_lpdf(y[i] | slope[Groups[i]]*x[i], sigma[Groups[i]]);
  }
}

