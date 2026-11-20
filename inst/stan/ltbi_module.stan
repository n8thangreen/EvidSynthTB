// Stage 1 LTBI Prevalence and Imperfect Diagnostic Test Module
// Supports:
// - Aggregate binomial survey data (M_agg studies/groups, with optional X_agg covariates)
// - Individual-level binary test data (N_ind individuals x K_tests tests)
// - Multiple diagnostic tests (latent class analysis when K_tests >= 2)
// - Optional covariates X_agg / X_ind on logit(pi)
// - Imperfect or fixed sensitivity and specificity

data {
  int<lower=1> K_tests;                       // Number of diagnostic tests
  int<lower=0> P_cov;                         // Number of covariates

  // Aggregate survey data
  int<lower=0> M_agg;
  array[M_agg] int<lower=0> y_agg;
  array[M_agg] int<lower=1> n_agg;
  array[M_agg] int<lower=1, upper=K_tests> test_id_agg;
  matrix[M_agg, P_cov] X_agg;

  // Individual-level diagnostic data
  int<lower=0> N_ind;
  matrix[N_ind, P_cov] X_ind;
  array[N_ind, K_tests] int<lower=0, upper=1> x_ind;

  // Imperfect test toggle and priors
  int<lower=0, upper=1> imperfect_test;
  vector<lower=0, upper=1>[K_tests] sens_fixed;
  vector<lower=0, upper=1>[K_tests] spec_fixed;
  vector<lower=0>[K_tests] a_sens;
  vector<lower=0>[K_tests] b_sens;
  vector<lower=0>[K_tests] a_spec;
  vector<lower=0>[K_tests] b_spec;

  // Priors on logit prevalence and regression coefficients
  real mu_alpha;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
}

parameters {
  real alpha_diag;
  vector[P_cov] beta_cov;
  vector<lower=0, upper=1>[imperfect_test ? K_tests : 0] sens_raw;
  vector<lower=0, upper=1>[imperfect_test ? K_tests : 0] spec_raw;
}

transformed parameters {
  real<lower=0, upper=1> prev_diag = inv_logit(alpha_diag);
  vector<lower=0, upper=1>[K_tests] sens;
  vector<lower=0, upper=1>[K_tests] spec;

  if (imperfect_test == 1) {
    sens = sens_raw;
    spec = spec_raw;
  } else {
    sens = sens_fixed;
    spec = spec_fixed;
  }
}

model {
  // Priors
  alpha_diag ~ normal(mu_alpha, sigma_alpha);
  if (P_cov > 0) {
    beta_cov ~ normal(0, sigma_beta);
  }
  if (imperfect_test == 1) {
    for (k in 1:K_tests) {
      sens_raw[k] ~ beta(a_sens[k], b_sens[k]);
      spec_raw[k] ~ beta(a_spec[k], b_spec[k]);
    }
  }

  // 1. Aggregate survey likelihood
  if (M_agg > 0) {
    vector[M_agg] eta_agg = rep_vector(alpha_diag, M_agg);
    if (P_cov > 0) {
      eta_agg += X_agg * beta_cov;
    }
    for (m in 1:M_agg) {
      int k = test_id_agg[m];
      real pi_m = inv_logit(eta_agg[m]);
      real p_pos = pi_m * sens[k] + (1.0 - pi_m) * (1.0 - spec[k]);
      y_agg[m] ~ binomial(n_agg[m], p_pos);
    }
  }

  // 2. Individual-level test likelihood (marginalizing over latent status Z_i)
  if (N_ind > 0) {
    vector[N_ind] eta = rep_vector(alpha_diag, N_ind);
    if (P_cov > 0) {
      eta += X_ind * beta_cov;
    }
    for (i in 1:N_ind) {
      real pi_i = inv_logit(eta[i]);
      real lp_z1 = log(pi_i);
      real lp_z0 = log1m(pi_i);
      for (k in 1:K_tests) {
        lp_z1 += bernoulli_lpmf(x_ind[i, k] | sens[k]);
        lp_z0 += bernoulli_lpmf(x_ind[i, k] | 1.0 - spec[k]);
      }
      target += log_sum_exp(lp_z0, lp_z1);
    }
  }
}

generated quantities {
  // Positive predictive value (PPV) for test 1 at baseline prevalence
  real<lower=0, upper=1> ppv_diag;
  {
    real p_pos1 = prev_diag * sens[1] + (1.0 - prev_diag) * (1.0 - spec[1]);
    ppv_diag = (prev_diag * sens[1]) / p_pos1;
  }
}
