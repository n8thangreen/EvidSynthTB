//
//
//

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

}

data {
  int<lower=0> N;
  vector[N] t;
  vector[N] d;

  //hyper parameters
  real mu_lambda;
  real sigma_lambda;
  real mu_gamma;
  real sigma_gamma;
}

parameters {
  real gamma;
  real<lower=0> lambda;
}

model {
  // priors
  lambda ~ normal(mu_lambda, sigma_lambda);
  gamma ~ normal(mu_gamma, sigma_gamma);

  // likelihood
  for (i in 1:N) {
    target += surv_gompertz_lpdf(t[i] | d[i], lambda, gamma);
  }
}

generated quantities {

  vector[10] ppred;
// what is j is not time
/// need t_pred[j]

  for (j in 1:10) {
    ppred[j] = gompertz_Surv(j, lambda, gamma);
  }

}
