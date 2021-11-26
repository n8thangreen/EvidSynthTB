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
}

data {
  int<lower=0> N;
  vector[N] t;
  vector[N] d;
  int<lower=0> t_lim;

  // hyper parameters
  real mu_lambda;
  real sigma_lambda;
  real mu_gamma;
  real sigma_gamma;

//TODO: alternative
  // int<lower=1> N_uncens;
  // int<lower=1> N_cens;
  // vector<lower=0>[N_cens] times_cens;
  // vector<lower=0>[N_uncens] times_uncens;
}

parameters {
  real loglambda;
  real loggamma;
}

transformed parameters {
  real<lower=0> lambda;
  real<lower=0> gamma;

  lambda = exp(loglambda);
  gamma = exp(loggamma);
}

model {

  // priors
  loglambda ~ normal(mu_lambda, sigma_lambda);
  loggamma ~ normal(mu_gamma, sigma_gamma);

  // likelihood
  for (i in 1:N) {
    target += surv_gompertz_lpdf(t[i] | d[i], gamma, lambda);
  }

//TODO: alternative formulation
    // target += gompertz_lpdf(times_uncens | gamma, lambda);
    // target += gompertz_lccdf(times_cens | gamma, lambda);
}

generated quantities {
  vector[t_lim] ppred;
// what is j is not time
/// need t_pred[j]

  for (j in 1:t_lim) {
    ppred[j] = gompertz_Surv(j, gamma, lambda);
  }

}

