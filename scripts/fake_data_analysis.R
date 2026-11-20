# Covariate-Adjusted Mixture Cure Model Analysis (Simulated Data)
# Simulates LTBI prevalence varying by standardized age (age_z) and high-incidence
# origin (high_inc_origin), and compares synthesis strategies with covariates.

library(EvidSynthTB)
library(dplyr)
library(ggplot2)

# 1. Simulate linked cohort with baseline covariates on logit(p_ltbi)
sim_cov <- create_fake_data(
  N_surv = 500,
  n_survey = 1000,
  p_ltbi = 0.30,
  sens = 0.88,
  spec = 0.96,
  surv_dist = "gompertz",
  shape = 0.15,
  rate = 0.10,
  include_covariates = TRUE,
  beta_cov = c(0.45, 0.70),
  seed = 401
)

# 2. Fit Joint (Latent-Class), Joint (Population), Cut, and Elastic models with covariates
comp_cov <- compare_evidsynth(
  prevalence_data = sim_cov$prevalence_data,
  progression_data = sim_cov$progression_data,
  models = c("joint", "cut", "elastic", "tempered"),
  coupling = "latent_class",
  imperfect_test = TRUE,
  covariates = c("age_z", "high_inc_origin"),
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 401
)

print(comp_cov$summary_table)
print(evaluate_comparison(comp_cov, truth = sim_cov$truth)$parameter_metrics)
print(diagnose_feedback(comp_cov))

# Inspect covariate regression coefficients beta_cov[1] and beta_cov[2]
posterior::summarise_draws(
  posterior::subset_draws(
    comp_cov$fits$joint$draws,
    variable = c("beta_cov[1]", "beta_cov[2]", "prev_cf", "lambda", "shape")
  )
)

plot_model_comparison(
  comp_cov,
  variables = c("prev_cf", "prev_diag", "lambda", "shape"),
  truth = c(prev_cf = 0.30, prev_diag = 0.30, lambda = 0.10, shape = 0.15)
)
plot_survival_comparison(comp_cov, type = "uncured")
