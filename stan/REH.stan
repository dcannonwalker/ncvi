//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N*G'.
data {
  int<lower=0> N; 
  int<lower=1> G;
  int<lower=1> P;
  int<lower=0> U;
  int<lower=0> y[N*G];
  int<lower=1, upper=G> group[N*G];
  row_vector[P] X[N*G];
  row_vector[U] Z[N*G];
  real<lower=0> sigB;
  real<lower=0> sigu;
  real<lower=0> sig;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu[P];
  vector[P] beta[G];
  vector[U] u[G];
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for (p in 1:P) {
    mu[p] ~ normal(0,sig);
    for (g in 1:G)
      beta[g,p] ~ normal(mu[p], sigB);
  }
  
  for (k in 1:U) {
    u[k] ~ normal(0, sigu);
  }
  
  for (n in 1:N*G)
    y[n] ~ poisson(exp(X[n]*beta[group[n]] + Z[n]*u[group[n]]));
}

