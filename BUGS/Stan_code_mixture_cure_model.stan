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
  int<lower=0> N_pos;      // ltbi positive
  int<lower=0> t_lim;

  // hyper parameters
  real mu_gamma;
  real sigma_gamma;
  real a_lambda;
  real b_lambda;

  // alternative
  int<lower=0> N_tb;
  int<lower=1> N_uncens;
  int<lower=1> N_cens;

  vector<lower=0>[N] t;
  vector<lower=0>[N_cens] t_cens;
  vector<lower=0>[N_uncens] t_uncens;

  // number of columns in the design matrix X
  int K;
  // design matrix X
  // should not include an intercept
  // matrix [N, K] X;

  vector[N] age;
  vector[N] ethnicity;

  // priors on regression coefficients
  real scale_alpha;
  vector[K] scale_beta;
}

parameters {
  real<lower=0> lambda;  // rate
  real gamma;            // shape

  // regression coefficient vector
  real alpha;
  vector[K] beta;
}

transformed parameters {
  vector[N] eta;
  real N_pos_hat;
  vector[N] p_pos;

  // eta = alpha + X * beta;
  eta = alpha + beta[1]*age + beta[2]*ethnicity;

  N_pos_hat = inv_logit(alpha)*N;
  p_pos = inv_logit(eta);
}

model {
  // priors
  lambda ~ normal(a_lambda, b_lambda);

  alpha ~ normal(0, scale_alpha);
  beta ~ normal(0, scale_beta);

  // likelihood

  // target = binomial_logit_lpmf(N_pos | N, eta);
  // for (k in 1:N) {
  //   // N_pos ~ bernoulli_logit(eta);
  //   target += bernoulli_logit_lpmf(pos[k] | eta[k]);
  // }


  // for (i in 1:N) {
  //   target += poisson_lpmf(N_tb | trisk*lambda);
  // }

  //// cure fraction model
  ///TODO: see bstanmcm code...

  for (i in 1:N) {
    target += log_sum_exp(
                log(cf[i]), log1m(cf[i]) + surv_gompertz_lpdf(t[i] | d[i], shape, lambda[i]));
  }

}

generated quantities {
  vector[t_lim] ppred;
  //TODO: what is j is not time
  //      need t_pred[j]

  for (j in 1:t_lim) {
    ppred[j] = gompertz_ccdf(j, gamma, lambda);
  }
}

