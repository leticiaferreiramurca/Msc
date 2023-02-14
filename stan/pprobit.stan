data {
  int<lower=0> k;
  int<lower=0> n;
  int<lower=0, upper=1> y[n];
  matrix[n,k] X;
}
parameters {
  real beta0;
  vector[k] beta;
  real delta;
}
transformed parameters {
  vector[n] prob;
  real lambda = exp(delta);

  for (i in 1:n) {
    real eta = (beta0 + X[i]*beta);
    prob[i] = pow(Phi(eta), lambda);
  }
}
model {
  beta0 ~ normal(0, 1000);
  beta ~ normal(0,1000);
  delta ~ uniform(-2,2);
  y ~ bernoulli(prob);
}

generated quantities {
  vector[n] log_lik;
  for (i in 1:n) {
    log_lik[i] = bernoulli_lpmf(y[i] | prob[i]);
  }
}
