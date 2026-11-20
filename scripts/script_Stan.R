# Reconstructed Abubakar (2018) UK PREDICT-TB Evidence Synthesis Analysis
# Combines baseline IGRA test positivity (1444 / 6380) with the digitized
# Kaplan-Meier individual participant survival data (raw-data/IPDdata.txt)
# across synthesis strategies and parametric survival families.

library(EvidSynthTB)
library(dplyr)
library(ggplot2)

abub <- load_abubakar_2018()

# 1. Compare synthesis strategies with imperfect IGRA sensitivity and specificity
comp_abub <- compare_evidsynth(
  prevalence_data = abub$prevalence_data,
  progression_data = abub$progression_data,
  models = c("cut", "joint", "independent", "melding", "elastic", "tempered"),
  imperfect_test = TRUE,
  surv_dist = "gompertz",
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 501
)

cat("=== Abubakar (2018) Synthesis Model Summary ===\n")
print(comp_abub$summary_table)

cat("=== Abubakar (2018) Inter-Module Feedback & Conflict Diagnostics ===\n")
print(diagnose_feedback(comp_abub))

# 2. Save fitted output for downstream ETS projection scripts
fit_cut <- comp_abub$fits$cut
s_cols <- grep("^S_pred\\[", names(fit_cut$draws), value = TRUE)
stan_output <- list(
  ppred = as.matrix(as.data.frame(fit_cut$draws)[, s_cols, drop = FALSE]),
  lambda = as.numeric(fit_cut$draws$lambda),
  gamma = as.numeric(fit_cut$draws$shape),
  prev_cf = as.numeric(fit_cut$draws$prev_cf)
)
if (dir.exists("data output")) {
  save(stan_output, comp_abub, file = file.path("data output", "stan_output.RData"))
}

# 3. Plots
plot_model_comparison(comp_abub, variables = c("prev_diag", "prev_cf", "lambda", "shape"))
plot_survival_comparison(comp_abub, type = "uncured")
plot_survival_comparison(comp_abub, type = "population")
