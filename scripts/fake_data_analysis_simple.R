# Baseline Simple Mixture Cure Model Analysis (Simulated Data)
# Compares all 6 synthesis strategies (joint, cut, independent, melding, elastic,
# tempered) and parametric survival families (gompertz, weibull, exponential)
# on simulated intercept-only LTBI and TB progression data.

library(EvidSynthTB)
library(dplyr)
library(ggplot2)

# 1. Simulate baseline artificial LTBI survey + progression cohort (no conflict)
sim <- create_fake_data(
  N_surv = 500,
  n_survey = 1000,
  p_ltbi = 0.30,
  sens = 1.0,
  spec = 1.0,
  surv_dist = "gompertz",
  shape = 0.1,
  rate = 0.1,
  t_censor = 10,
  seed = 101
)

# 2. Fit and compare all 6 synthesis strategies
comp_models <- compare_evidsynth(
  prevalence_data = sim$prevalence_data,
  progression_data = sim$progression_data[, c("t", "d")],
  models = c("joint", "cut", "independent", "melding", "elastic", "tempered"),
  imperfect_test = FALSE,
  surv_dist = "gompertz",
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 101
)

# 3. Quantitative comparison: recovery metrics & feedback diagnostics
eval_res <- evaluate_comparison(comp_models, truth = sim$truth)
print(eval_res$parameter_metrics)
print(eval_res$survival_metrics)

diag_res <- diagnose_feedback(comp_models)
print(diag_res)

# 4. Visual comparisons
plot_model_comparison(
  comp_models,
  variables = c("prev_cf", "prev_diag", "lambda", "shape"),
  truth = c(prev_cf = 0.30, prev_diag = 0.30, lambda = 0.10, shape = 0.10)
)
plot_survival_comparison(comp_models, type = "uncured")
plot_survival_comparison(comp_models, type = "population")

# 5. Compare parametric survival distributions (Gompertz vs Weibull vs Exponential)
fit_gomp <- comp_models$fits$joint
fit_weib <- evidsynth_fit(
  prevalence_data = sim$prevalence_data,
  progression_data = sim$progression_data[, c("t", "d")],
  model = "joint",
  surv_dist = "weibull",
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 101
)
fit_exp <- evidsynth_fit(
  prevalence_data = sim$prevalence_data,
  progression_data = sim$progression_data[, c("t", "d")],
  model = "joint",
  surv_dist = "exponential",
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 101
)

comp_dist <- structure(
  list(
    fits = list(gompertz = fit_gomp, weibull = fit_weib, exponential = fit_exp)
  ),
  class = "evidsynth_comparison"
)
print(diagnose_feedback(comp_dist))
plot_survival_comparison(comp_dist, type = "uncured")
