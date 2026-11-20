# Published Subgroup Summary Table Evidence Synthesis Analysis
# Replaces legacy JAGS/BUGS group-level count analysis with CmdStan mixture
# cure synthesis on published contingency tables (load_predict_aggregated()
# and raw-data/PREDICT_dummy.csv) without requiring individual-level data.

library(EvidSynthTB)
library(dplyr)
library(ggplot2)

# 1. Published multi-stratum PREDICT-TB contingency table (aggregated_data.csv)
pred_agg <- load_predict_aggregated(
  collapse_by = c("reasonforscreening", "inc_cob_participant2", "age_grp"),
  t_horizon = 5
)

print(head(pred_agg$summary_table))

# Fit and compare Cut, Joint, Elastic, and Tempered models directly on the
# published subgroup count table (n = pop, y = ltbi, d = tb) with covariates
comp_agg <- compare_evidsynth(
  prevalence_data = pred_agg$summary_table,
  progression_data = pred_agg$summary_table,
  models = c("cut", "joint", "elastic", "tempered"),
  imperfect_test = TRUE,
  covariates = c("reasonforscreening", "inc_cob_participant2"),
  surv_dist = "exponential",
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 801
)

cat("=== Published Subgroup Table Synthesis Summary ===\n")
print(comp_agg$summary_table)
print(diagnose_feedback(comp_agg))

plot_model_comparison(
  comp_agg,
  variables = c("prev_diag", "prev_cf", "lambda", "sens[1]", "spec[1]")
)
