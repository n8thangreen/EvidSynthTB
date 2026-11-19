test_that("load_abubakar_2018 and load_predict_migrant load real datasets cleanly", {
  abub <- load_abubakar_2018()
  expect_type(abub, "list")
  expect_named(abub, c("progression_data", "prevalence_data"))
  expect_equal(abub$prevalence_data$y, 1444L)
  expect_equal(abub$prevalence_data$n, 6380L)
  expect_gt(nrow(abub$progression_data), 6000)

  migrant <- load_predict_migrant()
  expect_type(migrant, "list")
  expect_named(migrant, c("progression_data", "prevalence_data"))
  expect_contains(
    names(migrant$progression_data),
    c("t", "d", "x", "age", "age_z", "ethnicity")
  )
  expect_gt(nrow(migrant$progression_data), 1000)
})
