# UKHSA / ETS Annual Surveillance Interval & Published KM Table Analysis
# Combines published LTBI survey tables with:
#   1. Annual active TB counts by year since UK entry (UKHSA / ETS style report)
#   2. Published Kaplan-Meier interval life tables (load_abubakar_km())
# across Cut, Joint, Independent, Melding, Elastic, and Tempered synthesis models.

library(EvidSynthTB)
library(dplyr)
library(ggplot2)

# 1. Published Kaplan-Meier Interval Table from Abubakar et al. (2018)
abub_km <- load_abubakar_km()
print(head(abub_km$progression_data))

comp_km <- compare_evidsynth(
  prevalence_data = abub_km$prevalence_data,
  progression_data = abub_km$progression_data,
  models = c("cut", "joint", "independent", "melding", "elastic", "tempered"),
  imperfect_test = TRUE,
  surv_dist = "gompertz",
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 901
)

cat("=== Abubakar (2018) Published KM Interval Table Synthesis ===\n")
print(comp_km$summary_table)
print(diagnose_feedback(comp_km))

# 2. UKHSA / ETS Annual Surveillance Report Synthesis (Year 1 to Year 5 post-entry)
pub_sim <- create_fake_published_data(
  n_groups = 6,
  n_per_group = 600,
  n_surveillance = 5000,
  p_ltbi = 0.28,
  p_ltbi_surv = 0.28,
  sens = 0.85,
  spec = 0.95,
  surv_dist = "gompertz",
  shape = -0.35,
  rate = 0.05,
  t_horizon = 5,
  seed = 902
)

print(pub_sim$interval_table)

comp_ukhsa <- compare_evidsynth(
  prevalence_data = pub_sim$prevalence_data,
  progression_data = pub_sim$interval_table,
  models = c("cut", "joint", "elastic", "tempered"),
  imperfect_test = TRUE,
  surv_dist = "gompertz",
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 902
)

cat("=== UKHSA / ETS Annual Interval Table Recovery ===\n")
print(evaluate_comparison(comp_ukhsa, truth = pub_sim$truth)$parameter_metrics)
print(diagnose_feedback(comp_ukhsa))

plot_survival_comparison(comp_ukhsa, type = "uncured")
