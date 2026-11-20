# Modular / Cut / Markov Melding Analysis on Real PREDICT-TB Migrant Cohort
# Compares Cut, Melding, Elastic, Tempered, Joint, and Independent synthesis
# models on the cleaned PREDICT-TB migrant dataset (load_predict_migrant()).

library(EvidSynthTB)
library(dplyr)
library(ggplot2)

migrant <- load_predict_migrant()

comp_migrant_modular <- compare_evidsynth(
  prevalence_data = migrant$prevalence_data,
  progression_data = migrant$progression_data,
  models = c("cut", "melding", "elastic", "tempered", "joint", "independent"),
  coupling = "population",
  imperfect_test = TRUE,
  covariates = "age_z",
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 701
)

cat("=== PREDICT-TB Migrant Cohort: Modular Synthesis Comparison ===\n")
print(comp_migrant_modular$summary_table)

cat("=== PREDICT-TB Migrant Cohort: Feedback & Conflict Diagnostics ===\n")
print(diagnose_feedback(comp_migrant_modular))

if (dir.exists("data output")) {
  save(
    comp_migrant_modular,
    file = file.path("data output", "stan_output_markov_melding.RData")
  )
}

plot_model_comparison(
  comp_migrant_modular,
  variables = c("prev_diag", "prev_cf", "lambda", "shape")
)
plot_survival_comparison(comp_migrant_modular, type = "uncured")
