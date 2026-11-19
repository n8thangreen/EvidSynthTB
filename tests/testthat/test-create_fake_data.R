test_that("create_fake_data generates well-formed synthetic datasets", {
  sim <- create_fake_data(
    N_surv = 200,
    n_survey = c(500, 800),
    p_ltbi = 0.30,
    p_ltbi_surv = 0.55,
    sens = c(0.85, 0.78),
    spec = c(0.95, 0.90),
    surv_dist = "gompertz",
    shape = 0.1,
    rate = 0.1,
    include_covariates = TRUE,
    seed = 42
  )

  expect_type(sim, "list")
  expect_named(sim, c("progression_data", "prevalence_data", "truth"))
  expect_equal(nrow(sim$progression_data), 200)
  expect_equal(nrow(sim$prevalence_data), 2)
  expect_contains(
    names(sim$progression_data),
    c("t", "d", "x", "x1", "x2", "z_true", "age_z", "high_inc_origin")
  )
  expect_equal(sim$truth$p_ltbi, 0.30)
  expect_equal(sim$truth$p_ltbi_surv, 0.55)
})

test_that("create_fake_data supports pos_only cohort and Weibull/Exponential families", {
  sim_pos <- create_fake_data(
    N_surv = 300,
    cohort_type = "pos_only",
    surv_dist = "weibull",
    shape = 1.2,
    rate = 0.08,
    seed = 99
  )

  expect_all_true(sim_pos$progression_data$x1 == 1L)
  expect_lt(nrow(sim_pos$progression_data), 300)
  expect_equal(sim_pos$truth$surv_dist, "weibull")
})
