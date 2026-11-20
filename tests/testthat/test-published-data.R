test_that("evidsynth_fit supports published summary tables, UKHSA interval tables, and KM life tables", {
  # 1. Simulated published summary table & UKHSA annual interval table
  pub_sim <- create_fake_published_data(
    n_groups = 4,
    n_per_group = 400,
    n_surveillance = 2000,
    p_ltbi = 0.25,
    p_ltbi_surv = 0.55,
    sens = 0.85,
    spec = 0.95,
    surv_dist = "gompertz",
    shape = -0.25,
    rate = 0.06,
    t_horizon = 5,
    seed = 42
  )

  expect_named(
    pub_sim,
    c("summary_table", "interval_table", "prevalence_data", "truth")
  )

  # Fit directly to annual surveillance interval table (Cut vs Joint)
  comp_interval <- compare_evidsynth(
    prevalence_data = pub_sim$prevalence_data,
    progression_data = pub_sim$interval_table,
    models = c("cut", "joint", "elastic"),
    imperfect_test = TRUE,
    surv_dist = "gompertz",
    chains = 1,
    n_iter = 300,
    n_burnin = 150,
    n_cut_subsamples = 6,
    seed = 111
  )

  expect_s3_class(comp_interval, "evidsynth_comparison")
  fb_int <- diagnose_feedback(comp_interval)
  expect_all_false(is.na(fb_int$waic_surv))
  expect_equal(fb_int$feedback_shift[fb_int$model == "cut"], 0, tolerance = 1e-8)

  # 2. Published Kaplan-Meier life table (load_abubakar_km)
  abub_km <- load_abubakar_km()
  expect_named(abub_km, c("progression_data", "prevalence_data"))
  expect_contains(
    names(abub_km$progression_data),
    c("t_start", "t_end", "n_risk", "d")
  )

  fit_km <- evidsynth_fit(
    prevalence_data = abub_km$prevalence_data,
    progression_data = abub_km$progression_data,
    model = "joint",
    imperfect_test = TRUE,
    surv_dist = "gompertz",
    chains = 1,
    n_iter = 300,
    n_burnin = 150,
    seed = 222
  )
  expect_s3_class(fit_km, "evidsynth_fit")

  # 3. Published multi-stratum PREDICT contingency table (load_predict_aggregated)
  pred_agg <- load_predict_aggregated(
    collapse_by = c("reasonforscreening", "inc_cob_participant2")
  )
  expect_contains(names(pred_agg$summary_table), c("n", "y", "d", "t_end"))

  fit_agg <- evidsynth_fit(
    progression_data = pred_agg$summary_table,
    model = "joint",
    imperfect_test = TRUE,
    covariates = "reasonforscreening",
    surv_dist = "exponential",
    chains = 1,
    n_iter = 300,
    n_burnin = 150,
    seed = 333
  )
  expect_s3_class(fit_agg, "evidsynth_fit")
})
