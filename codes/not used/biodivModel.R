
// defining data
data {
  int<lower=1> N;
  real x[N]; 
  real y[N];
}

parameters {
  real intercept;
  real slope;
  real<lower=0> sigma; 
}

model{
  for (i in 1:N)
  {
    target += normal_lpdf(y[i] | slope*x[i], sigma);
  }
  // Priors
  intercept ~ normal(0,0.001);
  slope ~ normal(0,0.001);
  sigma ~ normal(0,1)
}

