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
  real<lower=0> sigma_shape;
  real a_lambda;
  real b_lambda;
  real a_eth;
  real b_eth;

  vector<lower=0>[N] t;
  vector<lower=0, upper=1>[N] d;

  int pos[N];

  vector[N] age;
  int ethnicity[N];

  // priors on regression coefficients
  real<lower=0> scale_alpha;
  real<lower=0> scale_beta;
}

parameters {
  real<lower=0> lambda;  // rate
  real shape;

  // regression coefficient vector
  real alpha;
  real beta_age;

  vector[n_eth] beta_eth;
  real<lower=0> sd_eth;
}

transformed parameters {
  vector[N] eta;
  vector<lower=0, upper=1>[N] cf;

  eta = alpha + beta_age*age + beta_eth[ethnicity];
  cf = inv_logit(eta);
}

model {
  // priors
  lambda ~ gamma(a_lambda, b_lambda);
  shape ~ normal(mu_shape, sigma_shape);
  beta_age ~ normal(0, scale_beta);
  alpha ~ normal(0, scale_alpha);
  sd_eth ~ gamma(a_eth, b_eth);

  for (j in 1:n_eth) {
    beta_eth[j] ~ normal(0, sd_eth);
  }

////TODO: at diagnosis model?

  // likelihood
  // mixture cure model
  for (i in 1:N) {
    // pos[i] ~ bernoulli_logit(eta[i])

    // target += bernoulli_logit_lupmf(pos[i] | eta[i]) +
    //           log_sum_exp(
    //             log1m(cf[i]), log(cf[i]) + surv_gompertz_lpdf(t[i] | d[i], shape, lambda));

    target += log_sum_exp(
                log((cf[i]/(1-cf[i]))^pos[i]),
                log(cf[i]^(pos[i]+1) / (1 - cf[i])^(pos[i]-1)) + surv_gompertz_lpdf(t[i] | d[i], shape, lambda));
  }

}

generated quantities {
  vector[t_lim] S_pred;
  //TODO: what if j is not time
  //      need t_pred[j]

  for (j in 1:t_lim) {
    // S_pred[j] =(1 - inv_logit(alpha)) + inv_logit(alpha)*gompertz_Surv(j, shape, lambda);
    S_pred[j] = gompertz_Surv(j, shape, lambda);
  }
}

