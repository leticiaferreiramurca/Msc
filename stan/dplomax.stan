data {
  int<lower=0> k;
  int<lower=0> n;
  int<lower=0, upper=1> y[n];
  matrix[n,k] X;
}
parameters {
  real beta0;
  vector[k] beta;
  real loglambda;
}
transformed parameters {
  vector[n] prob;

  for (i in 1:n) {
    real eta = beta0 + X[i]*beta;
    if(eta > 0){
      prob[i] =  pow(1 - (0.5/ (1 + eta)), exp(loglambda));
    }else{
      prob[i] = pow((0.5 / (1 - eta)), exp(loglambda));
    }
  }
}
model {
  loglambda ~ uniform(-2,2);
  beta0 ~ normal(0, 100);
  beta ~ normal(0,100);
  y ~ bernoulli(prob);
}

generated quantities {
  vector[n] log_lik;
  for (i in 1:n) {
    log_lik[i] = bernoulli_lpmf(y[i] | prob[i]);
  }
}
