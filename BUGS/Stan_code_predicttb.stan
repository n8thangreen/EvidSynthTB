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
}

parameters {
  real mu;
  real<lower=0> lambda;
}

model {
  lambda ~ normal();
  mu ~ normal();

  // likelihood
  for (i in 1:N) {
    target += surv_gompertz_lpdf(t[i] | d[i], lambda[i], mu[i]);
  }
}

generated_quantities {

  vector[] pppred;

  for (j in 1:10) {
    ppred = surv_gompertz(t, lambda, mu);
  }

}
