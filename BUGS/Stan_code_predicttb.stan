// using predict-tb data
// estimate tb progression
// from ltbi positive cases

//TODO: how to vectorise likelihood functions?

functions {
/**
* gompertz
*
* @param t time
* @param rate
* @return A real
*/

  // // log hazard
  real gompertz_log_h (real t, real shape, real rate) {
    return(log(rate) + (shape * t));
  }

  // // hazard
  // real gompertz_haz (real t, real shape, real scale) {
  //   return(scale*exp(shape*t));
  // }

  // gompertz log survival
  real gompertz_lccdf (real t, real shape, real rate) {
    return(-rate/shape * (exp(shape * t) - 1));
  }

  // gompertz log density
  real gompertz_lpdf (real t, real shape, real rate) {
    return(gompertz_log_h(t, shape, rate) + gompertz_lccdf(t | shape, rate));
  }

  // gompertz survival
  real gompertz_ccdf (real t, real shape, real rate) {
    return(exp(-rate/shape * (exp(shape * t) - 1)));
  }

  // // gompertz log survival distribution
  // real surv_gompertz_lpdf (real t, real d, real shape, real rate) {
  //   return(d * gompertz_log_h(t, shape, rate) + gompertz_lccdf(t | shape, rate));
  // }
}

data {
  int<lower=0> N;
  int<lower=0> t_lim;

  // hyper parameters
  // real mu_lambda;
  // real sigma_lambda;
  real mu_gamma;
  real sigma_gamma;
  real a_lambda;
  real b_lambda;

  real mu0;
  real sigma0;
  real mu_eth;
  real sigma_eth;
  real mu_age;
  real sigma_age;

  // alternative
  int<lower=1> N_uncens;
  int<lower=1> N_cens;
  vector<lower=0>[N_cens] t_cens;
  vector<lower=0>[N_uncens] t_uncens;

  vector[N] eth;
  vector[N] age;
  int<lower=0,upper=1> pos[N];
}

parameters {
  // real loglambda;
  // real loggamma;
  real<lower=0> lambda;  // rate
  real gamma;   // shape

  real beta0;
  real beta_age;
  real beta_eth;
}

transformed parameters {
  // real<lower=0> lambda;
  //
  // lambda = exp(loglambda);
}

model {
  // priors
  // loglambda ~ normal(mu_lambda, sigma_lambda);

  lambda ~ gamma(a_lambda, b_lambda);
  gamma ~ normal(mu_gamma, sigma_gamma);

  beta0 ~ normal(mu0, sigma0);
  beta_eth ~ normal(mu_eth, sigma_eth);
  beta_age ~ normal(mu_age, sigma_age);

  // likelihood

  for (i in 1:N_uncens) {
    target += gompertz_lpdf(t_uncens[i] | gamma, lambda);
  }

  for (j in 1:N_cens) {
    target += gompertz_lccdf(t_cens[j] | gamma, lambda);
  }

  pos ~ bernoulli_logit(beta0 + beta_eth * eth + beta_age * age);
}

generated quantities {
  vector[t_lim] ppred;
  // what is j is not time
  // need t_pred[j]

  vector[89] age_tilde = 1:89;
  vector[6] eth_tilde = 1:6;
  matrix[89, 6]<lower=0,upper=1> y_tilde;

  for (j in 1:t_lim) {
    ppred[j] = gompertz_ccdf(j, gamma, lambda);
  }

  for (n in 1:89)
    for (m in 1:6)
      y_tilde[n, m] ~ bernoulli_logit_rng(beta0 + beta_eth * eth_tilde[m] + beta_age * age_tilde[n]);
    }
  }

}

