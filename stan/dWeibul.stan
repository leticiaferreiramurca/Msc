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

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> k;
  int<lower=0> n;
  int<lower=0, upper=1> y[n];
  matrix[n,k] X;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.

parameters {
  real beta0;
  vector[k] beta;
  real loglambda;
  real logc;
}


transformed parameters {
  real c = exp(logc);
  real lambda = exp(loglambda);
  vector[n] prob;
  
  for (i in 1:n) {
    real eta = -(beta0 + X[i]*beta);
    if(eta > 0){
      prob[i] = (1 - 0.5 * exp(-fabs(eta)^c))^lambda;
    }else{
       prob[i] = (0.5*exp(-fabs(eta)^c))^lambda;
    }
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  logc ~ normal(0,10);
  beta0 ~ normal(0.0,100);
  beta ~ normal(0.0,100);
  y ~ bernoulli(prob);
}

generated quantities {
  vector[n] log_lik;
  for (i in 1:n) {
    log_lik[i] = bernoulli_lpmf(y[i] | prob[i]);
  }
}
