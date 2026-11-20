# Individual-Level Linked PREDICT-TB Migrant Cohort Mixture Cure Analysis
# Uses load_predict_migrant() (cleaned_migrant_predict_data.RData) to fit
# covariate-adjusted mixture cure models comparing individual-level latent-class
# coupling vs population-level modular coupling.

library(EvidSynthTB)
library(dplyr)
library(ggplot2)

migrant <- load_predict_migrant()

# 1. Fit Joint model with individual-level latent-class coupling and age covariate
fit_lc <- evidsynth_fit(
  prevalence_data = migrant$prevalence_data,
  progression_data = migrant$progression_data,
  model = "joint",
  coupling = "latent_class",
  imperfect_test = TRUE,
  covariates = "age_z",
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 601
)

# 2. Fit Joint model with population-level modular coupling and age covariate
fit_pop <- evidsynth_fit(
  prevalence_data = migrant$prevalence_data,
  progression_data = migrant$progression_data,
  model = "joint",
  coupling = "population",
  imperfect_test = TRUE,
  covariates = "age_z",
  chains = 2,
  n_iter = 800,
  n_burnin = 300,
  seed = 601
)

comp_migrant_coupling <- structure(
  list(
    fits = list(latent_class = fit_lc, population = fit_pop),
    summary_table = bind_rows(
      mutate(fit_lc$summary, model = "latent_class"),
      mutate(fit_pop$summary, model = "population")
    ),
    draws_long = bind_rows(
      tibble::tibble(
        model = "latent_class",
        variable = "prev_cf",
        value = as.numeric(fit_lc$draws$prev_cf)
      ),
      tibble::tibble(
        model = "population",
        variable = "prev_cf",
        value = as.numeric(fit_pop$draws$prev_cf)
      )
    )
  ),
  class = "evidsynth_comparison"
)

print(comp_migrant_coupling$summary_table)
print(diagnose_feedback(comp_migrant_coupling))

# 3. Save output for downstream ETS migrant projections
s_cols <- grep("^S_pred\\[", names(fit_lc$draws), value = TRUE)
stan_output <- list(
  S_pred = as.matrix(as.data.frame(fit_lc$draws)[, s_cols, drop = FALSE]),
  ppred = as.matrix(as.data.frame(fit_lc$draws)[, s_cols, drop = FALSE]),
  lambda = as.numeric(fit_lc$draws$lambda),
  shape = as.numeric(fit_lc$draws$shape),
  cf = as.numeric(fit_lc$draws$prev_cf)
)
if (dir.exists("data output")) {
  save(
    stan_output,
    fit_lc,
    fit_pop,
    file = file.path("data output", "stan_output_mcm.RData")
  )
}

plot_progression(fit_lc, type = "uncured")
plot_survival_comparison(comp_migrant_coupling, type = "uncured")
