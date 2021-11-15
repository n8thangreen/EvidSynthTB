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

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  y ~ normal(mu, sigma);
  
    // likelihood
  for (i in 1:n_os) {
    target += log_sum_exp(
                log(cf_os) +
                surv_exp_lpdf(t_os[i] | d_os[i], lambda_os_bg[i]),
  }
}

