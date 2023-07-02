// estimate tb progression
// using kernel estimate posterior of latent prevalence
// modular models with feedback
// assume a perfect diagnostic test

// see:
// Goudie RJB, Presanis AM, Lunn D, De Angelis D, Wernisch L.
// Joining and splitting models with Markov melding. Bayesian Anal. 2019;14(1)
// p.95


functions {
#include /include/distributions.stan
}

data {
  int<lower=1> N;          // survival data sample size
  int<lower=0> t_lim;      // maximum time
  vector<lower=0>[N] t;
  vector<lower=0, upper=1>[N] d;  // censoring status

  // hyper parameters
  // gompertz
  real mu_shape;
  real<lower=0> sigma_shape;

  real a_lambda;  // rate
  real b_lambda;

  // cure fraction
  real mu_cf;
  real<lower=0> mu_sigma;

  // plug-in cure fraction data
  int<lower=1> M;          // diagnostic sample size
  real[M] mu_hat;
  real<lower=0>[M] sigma_hat;
}

parameters {
  real lin_cf;           // cure fraction linear equation
  real<lower=0> lambda;  // rate
  real shape;
}

transformed parameters {
  real<lower=0, upper=1> prev_cf;  // ltbi prevalance

  prev_cf = inv_logit(lin_cf);
}

model {
  // priors
  lambda ~ gamma(a_lambda, b_lambda);
  shape ~ normal(mu_shape, sigma_shape);

  lin_cf ~ normal(mu_cf, sigma_cf);

  // likelihood
  // mixture cure model (tb progression)
  for (i in 1:N) {
    target += log_sum_exp(
                log1m(prev_cf),
                log(prev_cf) + surv_gompertz_lpdf(t[i] | d[i], shape, lambda));
  }

  // incidence model (ltbi prevalence)
  for (i in 1:M) {
    target += normal_lpdf(mu_hat[i] | lin_cf, sigma_hat[i]);
  }
}

generated quantities {
  //TODO:
  // vector[t_lim] S_pred;
  //
  // for (j in 1:t_lim) {
  //   S_pred[j] = gompertz_Surv(j, shape, lambda);
  // }
}

