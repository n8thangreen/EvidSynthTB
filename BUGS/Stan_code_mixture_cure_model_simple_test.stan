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

  // hyper parameters
  real mu_shape;
  real<lower=0> sigma_shape;
  real a_lambda;
  real b_lambda;

  vector<lower=0>[N] t;
  vector<lower=0, upper=1>[N] d;

  int pos[N];
}

parameters {
  real<lower=0> lambda;  // rate
  real shape;
  real<lower=0, upper=1> sens;
  real<lower=0, upper=1> spec;
  real<lower=0, upper=1> prev_cf;
  real<lower=0, upper=1> prev_diag;
}

transformed parameters {
  real<lower=0, upper=1> prob_pos;
  prob_pos = prev_diag*sens + (1-prev_diag)*(1-spec);
}

model {
  // priors
  lambda ~ gamma(a_lambda, b_lambda);
  shape ~ normal(mu_shape, sigma_shape);
  sens ~ beta(100, 20);
  spec ~ beta(100, 20);
  prev_mean ~ normal(-1, 1);
  prev_sd ~ gamma(0.1, 0.1);

  # random effect
  logit(prev_diag) ~ normal(prev_mean, prev_sd);
  logit(prev_cf) ~ normal(prev_mean, prev_sd);

  // likelihood
  // mixture cure model
  for (i in 1:N) {
    pos[i] ~ bernoulli(prob_pos);

    target += log_sum_exp(
                log1m(prev_cf),
                log(prev_cf) + surv_gompertz_lpdf(t[i] | d[i], shape, lambda));
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

