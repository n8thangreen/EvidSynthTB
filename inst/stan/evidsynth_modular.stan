// Unified Bayesian Multi-Parameter Evidence Synthesis Model for EvidSynthTB
// Supports both:
//   - Published / Aggregate / Surveillance Report Data:
//     * Cumulative subgroup summary tables (N_g, y_g, d_g over [0, T_g])
//     * UKHSA / ETS annual / interval progression counts (d_g in [t0_g, t1_g] out of cohort N_g)
//     * Published Kaplan-Meier life-table intervals (d_g in [t0_g, t1_g] out of risk set R_g)
//   - Individual-Level Survival & Linked Diagnostic Data (t_i, d_i, x_i)
// And all 6 synthesis strategies:
//   1 = Joint, 2 = Independent, 3 = Cut Stage 2, 4 = Melding, 5 = Elastic, 6 = Tempered

functions {
#include distributions.stan
}

data {
  // Synthesis configuration
  int<lower=1, upper=6> model_type;
  int<lower=1, upper=2> coupling;
  int<lower=1, upper=2> cohort_type;
  int<lower=1, upper=3> surv_dist;
  int<lower=0, upper=1> imperfect_test;
  int<lower=1> t_lim;
  int<lower=0> P_cov;

  // Individual-level survival cohort data (Module 2a; N_surv = 0 if aggregate only)
  int<lower=0> N_surv;
  vector<lower=0>[N_surv] t;
  vector<lower=0, upper=1>[N_surv] d;
  matrix[N_surv, P_cov] X_surv;

  // Aggregate / interval / KM progression data (Module 2b; G_surv = 0 if IPD only)
  int<lower=0> G_surv;
  vector<lower=0>[G_surv] t0_agg;
  vector<lower=0>[G_surv] t1_agg;
  array[G_surv] int<lower=0> d_agg;
  array[G_surv] int<lower=1> n_risk_agg;
  array[G_surv] int<lower=0, upper=1> is_conditional_risk;
  matrix[G_surv, P_cov] X_surv_agg;

  // Diagnostic tests configuration
  int<lower=1> K_tests;
  int<lower=0, upper=1> is_linked;
  array[N_surv, K_tests] int<lower=0, upper=1> x_linked;

  // Aggregate survey LTBI testing data (Module 1a)
  int<lower=0> M_agg;
  array[M_agg] int<lower=0> y_agg;
  array[M_agg] int<lower=1> n_agg;
  array[M_agg] int<lower=1, upper=K_tests> test_id_agg;
  matrix[M_agg, P_cov] X_agg;

  // Unlinked individual-level diagnostic data (Module 1b)
  int<lower=0> N_ind;
  matrix[N_ind, P_cov] X_ind;
  array[N_ind, K_tests] int<lower=0, upper=1> x_ind;

  // Normal plug-in / Melding data (Module 1c)
  int<lower=0> M_plugin;
  vector[M_plugin] mu_hat;
  vector<lower=0>[M_plugin] sigma_hat;

  // Stage 1 fixed parameter values (used only when model_type == 3: Cut Stage 2)
  real fixed_alpha_diag;
  vector[P_cov] fixed_beta_cov;
  vector<lower=0, upper=1>[K_tests] fixed_sens;
  vector<lower=0, upper=1>[K_tests] fixed_spec;

  // Priors for diagnostic test sensitivity & specificity
  vector<lower=0, upper=1>[K_tests] sens_fixed;
  vector<lower=0, upper=1>[K_tests] spec_fixed;
  vector<lower=0>[K_tests] a_sens;
  vector<lower=0>[K_tests] b_sens;
  vector<lower=0>[K_tests] a_spec;
  vector<lower=0>[K_tests] b_spec;

  // Priors for prevalence and covariates
  real mu_alpha;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;

  // Priors for survival distribution (lambda = rate > 0, shape)
  real mu_shape;
  real<lower=0> sigma_shape;
  real<lower=0> a_lambda;
  real<lower=0> b_lambda;

  // Modular coupling hyperparameters
  real<lower=0> tau_couple;                   // SD for elastic coupling (model_type == 5)
  real<lower=0, upper=1> w_temper;            // Power-likelihood weight (model_type == 6)
}

transformed data {
  int has_alpha_diag = (model_type != 3);
  int has_separate_alpha_surv = (model_type == 2 || model_type == 5);
  int est_sens_spec = (model_type != 3 && imperfect_test == 1);
  int est_beta = (model_type != 3 && P_cov > 0);
  int est_shape = (surv_dist != 3);           // Exponential has no shape parameter
  int N_ll_total = N_surv + G_surv;
}

parameters {
  array[has_alpha_diag] real alpha_diag_raw;
  array[has_separate_alpha_surv] real alpha_surv_raw;
  vector[est_beta ? P_cov : 0] beta_cov_raw;

  vector<lower=0, upper=1>[est_sens_spec ? K_tests : 0] sens_raw;
  vector<lower=0, upper=1>[est_sens_spec ? K_tests : 0] spec_raw;

  real<lower=0> lambda;
  array[est_shape] real shape_raw;
}

transformed parameters {
  real alpha_diag;
  real alpha_surv;
  vector[P_cov] beta_cov;
  vector<lower=0, upper=1>[K_tests] sens;
  vector<lower=0, upper=1>[K_tests] spec;
  real shape;

  // 1. Resolve alpha_diag and beta_cov
  if (model_type == 3) {
    alpha_diag = fixed_alpha_diag;
    beta_cov = fixed_beta_cov;
    sens = fixed_sens;
    spec = fixed_spec;
  } else {
    alpha_diag = alpha_diag_raw[1];
    if (P_cov > 0) {
      beta_cov = beta_cov_raw;
    } else {
      beta_cov = rep_vector(0.0, 0);
    }
    if (imperfect_test == 1) {
      sens = sens_raw;
      spec = spec_raw;
    } else {
      sens = sens_fixed;
      spec = spec_fixed;
    }
  }

  // 2. Resolve alpha_surv
  if (has_separate_alpha_surv == 1) {
    alpha_surv = alpha_surv_raw[1];
  } else {
    alpha_surv = alpha_diag;
  }

  // 3. Resolve shape parameter
  if (surv_dist == 3) {
    shape = 0.0;
  } else if (surv_dist == 2) {
    shape = exp(shape_raw[1]);
  } else {
    shape = shape_raw[1];
  }

  real<lower=0, upper=1> prev_diag = inv_logit(alpha_diag);
  real<lower=0, upper=1> prev_cf = inv_logit(alpha_surv);
}

model {
  // --- Priors ---
  lambda ~ gamma(a_lambda, b_lambda);
  if (est_shape == 1) {
    shape_raw[1] ~ normal(mu_shape, sigma_shape);
  }

  if (model_type != 3) {
    alpha_diag_raw[1] ~ normal(mu_alpha, sigma_alpha);
    if (P_cov > 0) {
      beta_cov_raw ~ normal(0, sigma_beta);
    }
    if (imperfect_test == 1) {
      for (k in 1:K_tests) {
        sens_raw[k] ~ beta(a_sens[k], b_sens[k]);
        spec_raw[k] ~ beta(a_spec[k], b_spec[k]);
      }
    }
  }

  if (model_type == 2) {
    alpha_surv_raw[1] ~ normal(mu_alpha, sigma_alpha);
  } else if (model_type == 5) {
    alpha_surv_raw[1] ~ normal(alpha_diag, tau_couple);
  }

  // --- Module 1: LTBI Diagnostic / Prevalence Likelihood (when model_type != 3) ---
  if (model_type != 3) {
    // 1a. Aggregate binomial survey counts (with optional group covariates X_agg)
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

    // 1b. Unlinked individual diagnostic test results
    if (N_ind > 0) {
      vector[N_ind] eta_ind = rep_vector(alpha_diag, N_ind);
      if (P_cov > 0) {
        eta_ind += X_ind * beta_cov;
      }
      for (i in 1:N_ind) {
        real pi_i = inv_logit(eta_ind[i]);
        real lp_z1 = log(pi_i);
        real lp_z0 = log1m(pi_i);
        for (k in 1:K_tests) {
          lp_z1 += bernoulli_lpmf(x_ind[i, k] | sens[k]);
          lp_z0 += bernoulli_lpmf(x_ind[i, k] | 1.0 - spec[k]);
        }
        target += log_sum_exp(lp_z0, lp_z1);
      }
    }

    // 1c. Normal plug-in / Melding summary data
    for (m in 1:M_plugin) {
      target += normal_lpdf(mu_hat[m] | alpha_diag, sigma_hat[m]);
    }
  }

  // --- Module 2A: Individual-Level Progression Mixture Cure Likelihood ---
  if (N_surv > 0) {
    real w = (model_type == 6) ? w_temper : 1.0;
    vector[N_surv] eta_surv = rep_vector(alpha_surv, N_surv);
    vector[N_surv] eta_diag_surv = rep_vector(alpha_diag, N_surv);
    if (P_cov > 0) {
      eta_surv += X_surv * beta_cov;
      eta_diag_surv += X_surv * beta_cov;
    }

    if (is_linked == 1 && coupling == 1 && model_type != 3) {
      for (i in 1:N_surv) {
        real pi_diag_i = inv_logit(eta_diag_surv[i]);
        real lp_z1 = log(pi_diag_i);
        real lp_z0 = log1m(pi_diag_i);
        for (k in 1:K_tests) {
          lp_z1 += bernoulli_lpmf(x_linked[i, k] | sens[k]);
          lp_z0 += bernoulli_lpmf(x_linked[i, k] | 1.0 - spec[k]);
        }
        target += log_sum_exp(lp_z0, lp_z1);
      }
    }

    for (i in 1:N_surv) {
      real pi_surv_i = inv_logit(eta_surv[i]);
      real log_f = dispatch_lpdf(t[i] | surv_dist, shape, lambda);
      real log_S = dispatch_log_S(t[i], surv_dist, shape, lambda);

      if (is_linked == 1 && coupling == 2 && cohort_type == 1) {
        real lp_x_z1 = 0.0;
        real lp_x_z0 = 0.0;
        for (k in 1:K_tests) {
          lp_x_z1 += bernoulli_lpmf(x_linked[i, k] | sens[k]);
          lp_x_z0 += bernoulli_lpmf(x_linked[i, k] | 1.0 - spec[k]);
        }
        if (d[i] == 1) {
          if (model_type == 3) {
            target += log(pi_surv_i) + log_f;
          } else {
            target += lp_x_z1 + w * (log(pi_surv_i) + log_f);
          }
        } else {
          if (model_type == 3) {
            real log_num = log(pi_surv_i) + lp_x_z1;
            real log_den = log_sum_exp(log1m(pi_surv_i) + lp_x_z0, log_num);
            real ppv_i = exp(log_num - log_den);
            target += log_sum_exp(log1m(ppv_i), log(ppv_i) + log_S);
          } else {
            target += w * log_sum_exp(
              log1m(pi_surv_i) + lp_x_z0,
              log(pi_surv_i) + lp_x_z1 + log_S
            );
          }
        }
      } else {
        real pi_eff;
        if (cohort_type == 2) {
          real p_pos_i = pi_surv_i * sens[1] + (1.0 - pi_surv_i) * (1.0 - spec[1]);
          pi_eff = (pi_surv_i * sens[1]) / p_pos_i;
        } else {
          pi_eff = pi_surv_i;
        }

        if (d[i] == 1) {
          target += w * (log(pi_eff) + log_f);
        } else {
          target += w * log_sum_exp(log1m(pi_eff), log(pi_eff) + log_S);
        }
      }
    }
  }

  // --- Module 2B: Published Aggregate / Interval / KM Table Likelihood ---
  if (G_surv > 0) {
    real w = (model_type == 6) ? w_temper : 1.0;
    vector[G_surv] eta_agg_surv = rep_vector(alpha_surv, G_surv);
    if (P_cov > 0) {
      eta_agg_surv += X_surv_agg * beta_cov;
    }

    for (g in 1:G_surv) {
      real pi_surv_g = inv_logit(eta_agg_surv[g]);
      real pi_eff;
      if (cohort_type == 2) {
        real p_pos_g = pi_surv_g * sens[1] + (1.0 - pi_surv_g) * (1.0 - spec[1]);
        pi_eff = (pi_surv_g * sens[1]) / p_pos_g;
      } else {
        pi_eff = pi_surv_g;
      }

      real S0 = (t0_agg[g] <= 0.0) ? 1.0 : dispatch_Surv(t0_agg[g], surv_dist, shape, lambda);
      real S1 = dispatch_Surv(t1_agg[g], surv_dist, shape, lambda);
      real delta_S = fmax(S0 - S1, 1e-12);
      real prob_event;

      if (is_conditional_risk[g] == 1) {
        // Denominator is conditional risk set R_g at t0_agg[g] (e.g., KM life table)
        real S_pop_0 = (1.0 - pi_eff) + pi_eff * S0;
        prob_event = (pi_eff * delta_S) / fmax(S_pop_0, 1e-12);
      } else {
        // Denominator is initial cohort size N_g (e.g., cumulative table or UKHSA annual counts)
        prob_event = pi_eff * delta_S;
      }
      prob_event = fmin(fmax(prob_event, 1e-12), 1.0 - 1e-12);
      target += w * binomial_lpmf(d_agg[g] | n_risk_agg[g], prob_event);
    }
  }
}

generated quantities {
  vector[t_lim] S_pred;                       // Uncured survival S_u(t)
  vector[t_lim] S_pop_pred;                   // Population mixture survival (1 - prev_cf) + prev_cf * S_u(t)
  vector[N_ll_total] log_lik_surv;            // Pointwise survival log-likelihood for WAIC/LOO
  real<lower=0, upper=1> ppv_cf;

  {
    real p_pos1 = prev_cf * sens[1] + (1.0 - prev_cf) * (1.0 - spec[1]);
    ppv_cf = (prev_cf * sens[1]) / p_pos1;

    if (N_surv > 0) {
      vector[N_surv] eta_surv_gq = rep_vector(alpha_surv, N_surv);
      if (P_cov > 0) {
        eta_surv_gq += X_surv * beta_cov;
      }
      for (i in 1:N_surv) {
        real pi_surv_i = inv_logit(eta_surv_gq[i]);
        real log_f = dispatch_lpdf(t[i] | surv_dist, shape, lambda);
        real log_S = dispatch_log_S(t[i], surv_dist, shape, lambda);
        real pi_eff;

        if (is_linked == 1 && coupling == 2 && cohort_type == 1) {
          real lp_x_z1 = 0.0;
          real lp_x_z0 = 0.0;
          for (k in 1:K_tests) {
            lp_x_z1 += bernoulli_lpmf(x_linked[i, k] | sens[k]);
            lp_x_z0 += bernoulli_lpmf(x_linked[i, k] | 1.0 - spec[k]);
          }
          real log_num = log(pi_surv_i) + lp_x_z1;
          real log_den = log_sum_exp(log1m(pi_surv_i) + lp_x_z0, log_num);
          pi_eff = exp(log_num - log_den);
        } else if (cohort_type == 2) {
          real p_pos_i = pi_surv_i * sens[1] + (1.0 - pi_surv_i) * (1.0 - spec[1]);
          pi_eff = (pi_surv_i * sens[1]) / p_pos_i;
        } else {
          pi_eff = pi_surv_i;
        }

        if (d[i] == 1) {
          log_lik_surv[i] = log(pi_eff) + log_f;
        } else {
          log_lik_surv[i] = log_sum_exp(log1m(pi_eff), log(pi_eff) + log_S);
        }
      }
    }

    if (G_surv > 0) {
      vector[G_surv] eta_agg_gq = rep_vector(alpha_surv, G_surv);
      if (P_cov > 0) {
        eta_agg_gq += X_surv_agg * beta_cov;
      }
      for (g in 1:G_surv) {
        real pi_surv_g = inv_logit(eta_agg_gq[g]);
        real pi_eff;
        if (cohort_type == 2) {
          real p_pos_g = pi_surv_g * sens[1] + (1.0 - pi_surv_g) * (1.0 - spec[1]);
          pi_eff = (pi_surv_g * sens[1]) / p_pos_g;
        } else {
          pi_eff = pi_surv_g;
        }
        real S0 = (t0_agg[g] <= 0.0) ? 1.0 : dispatch_Surv(t0_agg[g], surv_dist, shape, lambda);
        real S1 = dispatch_Surv(t1_agg[g], surv_dist, shape, lambda);
        real delta_S = fmax(S0 - S1, 1e-12);
        real prob_event;
        if (is_conditional_risk[g] == 1) {
          real S_pop_0 = (1.0 - pi_eff) + pi_eff * S0;
          prob_event = (pi_eff * delta_S) / fmax(S_pop_0, 1e-12);
        } else {
          prob_event = pi_eff * delta_S;
        }
        prob_event = fmin(fmax(prob_event, 1e-12), 1.0 - 1e-12);
        log_lik_surv[N_surv + g] = binomial_lpmf(d_agg[g] | n_risk_agg[g], prob_event);
      }
    }
  }

  for (j in 1:t_lim) {
    S_pred[j] = dispatch_Surv(j * 1.0, surv_dist, shape, lambda);
    S_pop_pred[j] = (1.0 - prev_cf) + prev_cf * S_pred[j];
  }
}
