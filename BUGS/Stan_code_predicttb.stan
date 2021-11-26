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
  // real gompertz_log_h (real t, real shape, real scale) {
  //   return(log(scale) + (shape * t));
  // }
  //
  // // hazard
  // real gompertz_haz (real t, real shape, real scale) {
  //   return(scale*exp(shape*t));
  // }
  //
  // // gompertz log survival
  // real gompertz_lccdf (real t, real shape, real scale) {
  //   return(-scale/shape * (exp(shape * t) - 1));
  // }
  //
  // // gompertz log density
  // real gompertz_lpdf (real t, real shape, real scale) {
  //   return(gompertz_log_h(t, shape, scale) + gompertz_lccdf(t | shape, scale));
  // }
  //
  // // gompertz survival
  // real gompertz_ccdf (real t, real shape, real scale) {
  //   return(exp(-scale/shape * (exp(shape * t) - 1)));
  // }
  //
  // // gompertz log survival distribution
  // real surv_gompertz_lpdf (real t, real d, real shape, real scale) {
  //   return(d * gompertz_log_h(t, shape, scale) + gompertz_lccdf(t | shape, scale));
  // }
}

data {
  int<lower=0> N;
  int<lower=0> t_lim;

  // vector[N] t;
  // vector[N] d;

  // hyper parameters
  real mu_lambda;
  real sigma_lambda;
  real mu_gamma;
  real sigma_gamma;

  // alternative
  int<lower=1> N_uncens;
  int<lower=1> N_cens;
  vector<lower=0>[N_cens] t_cens;
  vector<lower=0>[N_uncens] t_uncens;
}

parameters {
  real loglambda;
  real loggamma;
}

transformed parameters {
  // real<lower=0> lambda;
  // real<lower=0> gamma;
  //
  // lambda = exp(loglambda);
  // gamma = exp(loggamma);
}

model {
  // // priors
  // loglambda ~ normal(mu_lambda, sigma_lambda);
  // loggamma ~ normal(mu_gamma, sigma_gamma);
  //
  // // likelihood
  //
  // // for (i in 1:N) {
  // //   target += surv_gompertz_lpdf(t[i] | d[i], gamma, lambda);
  // // }
  //
  // // alternative formulation
  // for (i in 1:N_uncens) {
  //   target += gompertz_lpdf(t_uncens[i] | gamma, lambda);
  // }
  //
  // for (j in 1:N_cens) {
  //   target += gompertz_lccdf(t_cens[j] | gamma, lambda);
  // }
}

generated quantities {
  // vector[t_lim] ppred;
  // // what is j is not time
  // // need t_pred[j]
  //
  // for (j in 1:t_lim) {
  //   ppred[j] = gompertz_ccdf(j, gamma, lambda);
  // }

}
