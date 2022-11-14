// using predict-tb data
// estimate tb progression
// from ltbi positive cases
// hierarchical model on prevalence

functions {
#include /include/distributions.stan
}

data {
  int<lower=1> N;          // total sample size
  int<lower=1> N_obs;
  int<lower=1> N_cens;
  int<lower=0> t_lim;

  int obs_idx[N_obs];
  int cens_idx[N_cens];

  // hyper parameters
  real mu_shape;
  real<lower=0> sigma_shape;
  real a_lambda;
  real b_lambda;

  vector<lower=0>[N] t;
  vector<lower=0, upper=1>[N] d;

  int x[N];  // positive test results
}

parameters {
  real<lower=0, upper=1> sens;
  real<lower=0, upper=1> spec;
  real lin_diag;
  real lin_cf;
  real<lower=0> lambda;  // rate
  real shape;
  real prev_mean;
  real<lower=0> prev_sd;
}

transformed parameters {
  real<lower=0, upper=1> prob_pos;
  real<lower=0, upper=1> prev_diag;
  real<lower=0, upper=1> prev_cf;

  prev_diag = inv_logit(lin_diag);
  prev_cf = inv_logit(lin_cf);

  prob_pos = prev_diag*sens + (1-prev_diag)*(1-spec);
}

model {
  // priors
  lambda ~ gamma(a_lambda, b_lambda);
  shape ~ normal(mu_shape, sigma_shape);

  // likelihood
  // mixture cure model
  for (i in 1:N) {
    target += log_sum_exp(
                log1m(prev_cf),
                log(prev_cf) + surv_gompertz_lpdf(t[i] | d[i], shape, lambda));
  }

  for (i in 1:N_cens) {
    x[cens_idx[i]] ~ bernoulli(prob_pos);
  }

  for (i in 1:N_obs) {
    x[obs_idx[i]] ~ bernoulli(sens);
  }
}

generated quantities {
  vector[t_lim] S_pred;
  //TODO: what if j is not time
  //      need t_pred[j]

  for (j in 1:t_lim) {
    S_pred[j] = gompertz_Surv(j, shape, lambda);
  }
}

