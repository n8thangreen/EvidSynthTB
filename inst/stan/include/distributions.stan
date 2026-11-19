// Survival distribution helper functions for EvidSynthTB
// Compatible with CmdStan 2.33+

real normal_pdf(real t, real mu, real sigma) {
  return 1.0 / (sigma * sqrt(2.0 * pi())) * exp(-0.5 * pow((t - mu) / sigma, 2));
}

// ============================================================================
// 1. Exponential distribution
// ============================================================================

real exp_log_h(real t, real rate) {
  return log(rate);
}

real exp_haz(real t, real rate) {
  return rate;
}

real exp_log_S(real t, real rate) {
  return -rate * t;
}

real exp_Surv(real t, real rate) {
  return exp(-rate * t);
}

real surv_exp_pdf(real t, real d, real rate) {
  return pow(exp_haz(t, rate), d) * exp_Surv(t, rate);
}

real surv_exp_lpdf(real t, real d, real rate) {
  return d * exp_log_h(t, rate) + exp_log_S(t, rate);
}

// ============================================================================
// 2. Weibull distribution (shape, rate/scale parameterization)
//    Matches flexsurv::rweibull(shape, scale = 1 / rate) or scale directly
//    Here 'rate' acts as the hazard scale parameter lambda > 0:
//    H(t) = rate * t^shape, h(t) = rate * shape * t^(shape - 1)
// ============================================================================

real weibull_ph_log_h(real t, real shape, real rate) {
  return log(rate) + log(shape) + (shape - 1.0) * log(t);
}

real weibull_ph_log_S(real t, real shape, real rate) {
  return -rate * pow(t, shape);
}

real weibull_ph_Surv(real t, real shape, real rate) {
  return exp(-rate * pow(t, shape));
}

// Traditional AFT Weibull (shape, scale) preserved for compatibility
real weibull_log_h(real t, real shape, real scale) {
  return log(shape) + (shape - 1.0) * log(t / scale) - log(scale);
}

real weibull_haz(real t, real shape, real scale) {
  return shape / scale * pow(t / scale, shape - 1.0);
}

real weibull_log_S(real t, real shape, real scale) {
  return -pow(t / scale, shape);
}

real weibull_Surv(real t, real shape, real scale) {
  return exp(-pow(t / scale, shape));
}

real surv_weibull_lpdf(real t, real d, real shape, real scale) {
  return d * weibull_log_h(t, shape, scale) + weibull_log_S(t, shape, scale);
}

// ============================================================================
// 3. Gompertz distribution (matches flexsurv::pgompertz(shape, rate))
//    h(t) = rate * exp(shape * t)
//    S(t) = exp(-rate / shape * (exp(shape * t) - 1))
// ============================================================================

real gompertz_log_h(real t, real shape, real rate) {
  return log(rate) + (shape * t);
}

real gompertz_haz(real t, real shape, real rate) {
  return rate * exp(shape * t);
}

real gompertz_log_S(real t, real shape, real rate) {
  if (abs(shape) < 1e-8) {
    return -rate * t;
  }
  return -rate / shape * (exp(shape * t) - 1.0);
}

real gompertz_Surv(real t, real shape, real rate) {
  return exp(gompertz_log_S(t, shape, rate));
}

real gompertz_lpdf(real t, real shape, real rate) {
  return gompertz_log_h(t, shape, rate) + gompertz_log_S(t, shape, rate);
}

real surv_gompertz_lpdf(real t, real d, real shape, real rate) {
  return d * gompertz_log_h(t, shape, rate) + gompertz_log_S(t, shape, rate);
}

real inv_cdf_gompertz(real p, real shape, real rate) {
  return 1.0 / shape * log(1.0 + (-shape / rate * log1m(p)));
}

// ============================================================================
// 4. Unified survival dispatcher
//    dist_id: 1 = Gompertz, 2 = Weibull (PH), 3 = Exponential
// ============================================================================

real dispatch_log_h(real t, int dist_id, real shape, real rate) {
  if (dist_id == 1) {
    return gompertz_log_h(t, shape, rate);
  } else if (dist_id == 2) {
    return weibull_ph_log_h(t, shape, rate);
  } else {
    return exp_log_h(t, rate);
  }
}

real dispatch_log_S(real t, int dist_id, real shape, real rate) {
  if (dist_id == 1) {
    return gompertz_log_S(t, shape, rate);
  } else if (dist_id == 2) {
    return weibull_ph_log_S(t, shape, rate);
  } else {
    return exp_log_S(t, rate);
  }
}

real dispatch_lpdf(real t, int dist_id, real shape, real rate) {
  return dispatch_log_h(t, dist_id, shape, rate) + dispatch_log_S(t, dist_id, shape, rate);
}

real dispatch_Surv(real t, int dist_id, real shape, real rate) {
  return exp(dispatch_log_S(t, dist_id, shape, rate));
}

// ============================================================================
// 5. Log-logistic, Log-normal, and Generalized Gamma
// ============================================================================

real loglogistic_log_h(real t, real shape, real scale) {
  return log(shape) - log(scale) + (shape - 1.0) * (log(t) - log(scale)) - log1p(pow(t / scale, shape));
}

real loglogistic_haz(real t, real shape, real scale) {
  return (shape / scale * pow(t / scale, shape - 1.0)) / (1.0 + pow(t / scale, shape));
}

real loglogistic_log_S(real t, real shape, real scale) {
  return -log1p(pow(t / scale, shape));
}

real loglogistic_Surv(real t, real shape, real scale) {
  return 1.0 / (1.0 + pow(t / scale, shape));
}

real surv_loglogistic_lpdf(real t, real d, real shape, real scale) {
  return d * loglogistic_log_h(t, shape, scale) + loglogistic_log_S(t, shape, scale);
}

real gen_gamma_lpdf(real t, real mu, real sigma, real Q) {
  real w = (log(t) - mu) / sigma;
  real q_inv2 = pow(Q, -2);
  return -log(sigma * t) + log(abs(Q)) + q_inv2 * log(q_inv2) + q_inv2 * (Q * w - exp(Q * w)) - lgamma(q_inv2);
}

real gen_gamma_Surv(real t, real mu, real sigma, real Q) {
  real w = (log(t) - mu) / sigma;
  real qq = 1.0 / (Q * Q);
  real expnu = exp(Q * w) * qq;
  return 1.0 - gamma_cdf(expnu | qq, 1.0);
}

real gen_gamma_log_S(real t, real mu, real sigma, real Q) {
  return log(gen_gamma_Surv(t, mu, sigma, Q));
}

real gen_gamma_log_h(real t, real mu, real sigma, real Q) {
  return gen_gamma_lpdf(t | mu, sigma, Q) - gen_gamma_log_S(t, mu, sigma, Q);
}

real gen_gamma_haz(real t, real mu, real sigma, real Q) {
  return exp(gen_gamma_log_h(t, mu, sigma, Q));
}

real lognormal_log_S(real t, real mu, real sigma) {
  return log1m(Phi((log(t) - mu) / sigma));
}

real lognormal_log_h(real t, real mu, real sigma) {
  return lognormal_lpdf(t | mu, sigma) - lognormal_lccdf(t | mu, sigma);
}

real lognormal_Surv(real t, real mu, real sigma) {
  return 1.0 - Phi((log(t) - mu) / sigma);
}

real lognormal_haz(real t, real mu, real sigma) {
  return (1.0 / (t * sigma)) * normal_pdf(log(t), mu, sigma) / (1.0 - Phi((log(t) - mu) / sigma));
}

real surv_lognormal_lpdf(real t, real d, real mu, real sigma) {
  return d * lognormal_log_h(t, mu, sigma) + lognormal_log_S(t, mu, sigma);
}
