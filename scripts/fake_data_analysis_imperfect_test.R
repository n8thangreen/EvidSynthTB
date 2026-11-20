# Imperfect Diagnostic Test & Mixture Cure Model Analysis (Simulated Data)
# Evaluates:
#   1. Population-level coupling vs Individual-level latent-class coupling
#   2. Single imperfect test (K = 1) vs Dual imperfect tests (K = 2, e.g. IGRA + TST)
#   3. All-comers follow-up (cohort_type = "all") vs Test-positive subgroup
#      follow-up (cohort_type = "pos_only", where uncured fraction is PPV)

library(EvidSynthTB)
library(dplyr)
library(ggplot2)

# 1. Simulate linked cohort with K = 2 imperfect diagnostic tests (IGRA + TST)
sim_k2 <- create_fake_data(
  N_surv = 600,
  n_survey = 1200,
  p_ltbi = 0.30,
  sens = c(0.85, 0.78),
  spec = c(0.95, 0.90),
  surv_dist = "gompertz",
  shape = 0.2,
  rate = 0.10,
  t_censor = 10,
  seed = 301
)

# 2. Compare Population-level vs Individual-level Latent-Class Coupling (K = 1 vs K = 2)
fit_pop_k1 <- evidsynth_fit(
  prevalence_data = sim_k2$prevalence_data,
  progression_data = sim_k2$progression_data[, c("t", "d", "x1")],
  model = "joint",
  coupling = "population",
  imperfect_test = TRUE,
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 301
)

fit_lc_k1 <- evidsynth_fit(
  prevalence_data = sim_k2$prevalence_data,
  progression_data = sim_k2$progression_data[, c("t", "d", "x1")],
  model = "joint",
  coupling = "latent_class",
  imperfect_test = TRUE,
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 301
)

fit_lc_k2 <- evidsynth_fit(
  prevalence_data = sim_k2$prevalence_data,
  progression_data = sim_k2$progression_data[, c("t", "d", "x1", "x2")],
  model = "joint",
  coupling = "latent_class",
  imperfect_test = TRUE,
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 301
)

comp_coupling <- structure(
  list(
    fits = list(
      pop_1test = fit_pop_k1,
      latent_class_1test = fit_lc_k1,
      latent_class_2tests = fit_lc_k2
    ),
    summary_table = bind_rows(
      mutate(fit_pop_k1$summary, model = "pop_1test"),
      mutate(fit_lc_k1$summary, model = "latent_class_1test"),
      mutate(fit_lc_k2$summary, model = "latent_class_2tests")
    )
  ),
  class = "evidsynth_comparison"
)

cat("=== Coupling & Multi-Test Recovery ===\n")
print(evaluate_comparison(comp_coupling, truth = sim_k2$truth)$parameter_metrics)
print(diagnose_feedback(comp_coupling))

# 3. Test-Positive Subgroup Follow-Up (cohort_type = "pos_only", uncured fraction = PPV)
sim_pos <- create_fake_data(
  N_surv = 800,
  n_survey = 1500,
  p_ltbi = 0.30,
  sens = 0.85,
  spec = 0.92,
  surv_dist = "gompertz",
  shape = 0.15,
  rate = 0.10,
  cohort_type = "pos_only",
  seed = 302
)

comp_pos <- compare_evidsynth(
  prevalence_data = sim_pos$prevalence_data,
  progression_data = sim_pos$progression_data[, c("t", "d")],
  models = c("joint", "cut", "elastic"),
  cohort_type = "pos_only",
  imperfect_test = TRUE,
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 302
)

cat("=== Test-Positive Subgroup Cohort (PPV Cure Fraction) ===\n")
print(evaluate_comparison(comp_pos, truth = sim_pos$truth)$parameter_metrics)
plot_model_comparison(
  comp_pos,
  variables = c("prev_diag", "ppv_cf", "lambda", "shape"),
  truth = c(
    prev_diag = sim_pos$truth$p_ltbi,
    ppv_cf = sim_pos$truth$ppv_surv,
    lambda = sim_pos$truth$rate,
    shape = sim_pos$truth$shape
  )
)
