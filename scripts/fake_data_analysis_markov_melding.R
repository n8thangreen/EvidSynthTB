# Modular Synthesis, Cut Models, and Markov Melding Analysis (Simulated Data)
# Compares all 6 synthesis strategies under:
#   (A) Datasets in Agreement (p_ltbi = 0.30, p_ltbi_surv = 0.30)
#   (B) Datasets in Conflict / Selection Bias (p_ltbi = 0.25, p_ltbi_surv = 0.65)
# Demonstrates how Cut, Elastic, and Tempered models control feedback from a
# biased progression cohort into LTBI prevalence estimates.

library(EvidSynthTB)
library(dplyr)
library(ggplot2)

# --- Scenario A: Datasets in Agreement ---
sim_agree <- create_fake_data(
  N_surv = 500,
  n_survey = c(1000, 1000),
  p_ltbi = 0.30,
  p_ltbi_surv = 0.30,
  sens = 0.85,
  spec = 0.95,
  surv_dist = "gompertz",
  shape = 0.1,
  rate = 0.12,
  seed = 202
)

comp_agree <- compare_evidsynth(
  prevalence_data = sim_agree$prevalence_data,
  progression_data = sim_agree$progression_data[, c("t", "d")],
  models = c("joint", "cut", "independent", "melding", "elastic", "tempered"),
  imperfect_test = TRUE,
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 202
)

cat("=== Scenario A (Agreement): Feedback & Conflict Diagnostics ===\n")
print(diagnose_feedback(comp_agree))
print(evaluate_comparison(comp_agree, truth = sim_agree$truth)$parameter_metrics)

# --- Scenario B: Datasets in Conflict (Biased Progression Cohort) ---
sim_conflict <- create_fake_data(
  N_surv = 500,
  n_survey = c(1000, 1000),
  p_ltbi = 0.25,
  p_ltbi_surv = 0.65,
  sens = 0.85,
  spec = 0.95,
  surv_dist = "gompertz",
  shape = 0.1,
  rate = 0.12,
  seed = 203
)

comp_conflict <- compare_evidsynth(
  prevalence_data = sim_conflict$prevalence_data,
  progression_data = sim_conflict$progression_data[, c("t", "d")],
  models = c("joint", "cut", "independent", "melding", "elastic", "tempered"),
  imperfect_test = TRUE,
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 203
)

cat("=== Scenario B (Conflict): Feedback & Conflict Diagnostics ===\n")
print(diagnose_feedback(comp_conflict))
print(evaluate_comparison(comp_conflict, truth = sim_conflict$truth)$parameter_metrics)

# Save output and plot posterior comparisons
if (dir.exists("data output")) {
  save(
    comp_agree,
    comp_conflict,
    file = file.path("data output", "stan_output_fake_markov_melding.RData")
  )
}

plot_model_comparison(
  comp_conflict,
  variables = c("prev_diag", "prev_cf", "lambda", "shape"),
  truth = c(prev_diag = 0.25, prev_cf = 0.65, lambda = 0.12, shape = 0.10)
)
plot_survival_comparison(comp_conflict, type = "uncured")
