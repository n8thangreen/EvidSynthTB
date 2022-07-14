// using predict-tb data
// estimate tb progression
// from ltbi positive cases


functions {
/**
* gompertz
*
* @param t time
* @param rate
* @return A real
*/

// log hazard
real gompertz_log_h (real t, real shape, real scale) {
  real log_h;
  log_h = log(scale) + (shape * t);
  return log_h;
}

// hazard
real gompertz_haz (real t, real shape, real scale) {
  real h;
  h = scale*exp(shape*t);
  return h;
}

// gompertz log survival
real gompertz_log_S (real t, real shape, real scale) {
  real log_S;
  log_S = -scale/shape * (exp(shape * t) - 1);
  return log_S;
}

// gompertz survival
real gompertz_Surv (real t, real shape, real scale) {
  real S;
  S = exp(-scale/shape * (exp(shape * t) - 1));
  return S;
}

// gompertz sampling distribution
real surv_gompertz_lpdf (real t, real d, real shape, real scale) {
  real log_lik;
  log_lik = d * gompertz_log_h(t, shape, scale) + gompertz_log_S(t, shape, scale);
  return log_lik;
}

// quantile
real inv_cdf_gompertz (real p, real shape, real scale) {
  real res;
  res = 1/shape * log(1 + (-shape/scale * log(1 - p)));
  return res;
}
}

data {
  int<lower=1> N;          // total sample size
  int<lower=0> t_lim;
  int<lower=0> n_eth;

  // hyper parameters
  real mu_shape;
  real sigma_shape;
  real a_lambda;
  real b_lambda;
  real a_eth;
  real b_eth;

  vector<lower=0>[N] t;
  vector<lower=0>[N] d;

  vector[N] pos;

  vector[N] age;
  int ethnicity[N];

  // priors on regression coefficients
  real scale_alpha;
  real scale_beta;
}

parameters {
  real<lower=0> lambda;  // rate
  real shape;

  // regression coefficient vector
  real alpha;
  real beta_age;

  vector[n_eth] beta_eth;
  real sd_eth;
}

transformed parameters {
  vector[N] eta;
  vector[N] cf;

  eta = alpha + beta_age*age + beta_eth[ethnicity];
  cf = inv_logit(eta);
}

model {
  // priors
  lambda ~ normal(a_lambda, b_lambda);
  shape ~ normal(mu_shape, sigma_shape);
  beta_age ~ normal(0, scale_beta);
  alpha ~ normal(0, scale_alpha);
  sd_eth ~ gamma(a_eth, b_eth);

  for (j in 1:n_eth) {
    beta_eth[j] ~ normal(0, sd_eth);
  }

  for (i in 1:N) {
    pos[i] ~ bernoulli(cf[i])
  }

  // likelihood
  // cure fraction model
  for (i in 1:N) {
    target += log_sum_exp(
                log(cf[i]), log1m(cf[i]) + surv_gompertz_lpdf(t[i] | d[i], shape, lambda));
  }

}

generated quantities {
  vector[t_lim] S_pred;
  //TODO: what is j is not time
  //      need t_pred[j]

  for (j in 1:t_lim) {
    S_pred[j] = gompertz_Surv(j, shape, lambda);
  }
}

