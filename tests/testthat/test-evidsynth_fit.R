test_that("evidsynth_fit works across synthesis strategies and blocks feedback in cut model", {
  # Simulate conflict scenario: survey true prevalence = 0.25, survival cohort = 0.70
  sim_no_conflict <- create_fake_data(
    N_surv = 150,
    n_survey = 2000,
    p_ltbi = 0.25,
    p_ltbi_surv = 0.25,
    sens = 0.85,
    spec = 0.95,
    shape = 0.1,
    rate = 0.15,
    seed = 101
  )

  sim_conflict <- create_fake_data(
    N_surv = 150,
    n_survey = 2000,
    p_ltbi = 0.25,
    p_ltbi_surv = 0.75,
    sens = 0.85,
    spec = 0.95,
    shape = 0.1,
    rate = 0.15,
    seed = 102
  )
  # Use identical survey data in both to verify cut invariance
  sim_conflict$prevalence_data <- sim_no_conflict$prevalence_data

  # Drop individual x from progression_data so Module 1 is strictly the survey
  prog_clean <- sim_no_conflict$progression_data[, c("t", "d")]
  prog_biased <- sim_conflict$progression_data[, c("t", "d")]

  fit_cut_clean <- evidsynth_fit(
    prevalence_data = sim_no_conflict$prevalence_data,
    progression_data = prog_clean,
    model = "cut",
    imperfect_test = TRUE,
    chains = 1,
    n_iter = 300,
    n_burnin = 150,
    n_cut_subsamples = 8,
    seed = 555
  )

  fit_cut_biased <- evidsynth_fit(
    prevalence_data = sim_conflict$prevalence_data,
    progression_data = prog_biased,
    model = "cut",
    imperfect_test = TRUE,
    chains = 1,
    n_iter = 300,
    n_burnin = 150,
    n_cut_subsamples = 8,
    seed = 555
  )

  # In a true cut model, prev_diag and prev_cf draws must be identical regardless
  # of biased progression survival data!
  expect_equal(
    mean(fit_cut_clean$draws$prev_diag),
    mean(fit_cut_biased$draws$prev_diag),
    tolerance = 1e-6
  )

  # Compare models and plotting helpers
  comp <- compare_evidsynth(
    prevalence_data = sim_conflict$prevalence_data,
    progression_data = prog_biased,
    models = c("joint", "independent", "elastic", "tempered"),
    imperfect_test = TRUE,
    chains = 1,
    n_iter = 300,
    n_burnin = 150,
    seed = 777
  )

  expect_s3_class(comp, "evidsynth_comparison")
  expect_named(comp$fits, c("joint", "independent", "elastic", "tempered"))

  p1 <- plot_progression(comp$fits$joint)
  p2 <- plot_model_comparison(comp)
  p3 <- plot_survival_comparison(comp)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})

test_that("evidsynth_fit supports linked latent_class coupling with imperfect tests", {
  sim <- create_fake_data(
    N_surv = 120,
    n_survey = 500,
    p_ltbi = 0.35,
    sens = c(0.85, 0.80),
    spec = c(0.95, 0.92),
    include_covariates = TRUE,
    seed = 202
  )

  fit_lc <- evidsynth_fit(
    prevalence_data = sim$prevalence_data,
    progression_data = sim$progression_data,
    model = "joint",
    coupling = "latent_class",
    imperfect_test = TRUE,
    covariates = c("age_z", "high_inc_origin"),
    chains = 1,
    n_iter = 300,
    n_burnin = 150,
    seed = 303
  )

  expect_s3_class(fit_lc, "evidsynth_fit")
  expect_contains(
    fit_lc$summary$variable,
    c("prev_diag", "prev_cf", "lambda", "shape", "sens[1]", "spec[1]")
  )
})
