#' Load reconstructed Abubakar (2018) UK PREDICT-TB data
#'
#' Loads the digitized individual participant survival data (`IPDdata.txt`)
#' reconstructed from the Kaplan-Meier curve in Abubakar et al. (2018) alongside
#' the aggregate baseline IGRA positivity count (`y = 1444` positive out of
#' `n = 6380` tested).
#'
#' @return A list with `progression_data` (tibble with columns `t` and `d`) and
#'   `prevalence_data` (tibble with `y = 1444L`, `n = 6380L`, `test_id = 1L`,
#'   and empirical logit plug-in estimates `mu_hat` and `sigma_hat`).
#' @export
load_abubakar_2018 <- function() {
  ipd_path <- system.file("extdata", "IPDdata.txt", package = "EvidSynthTB")
  if (!nzchar(ipd_path) || !file.exists(ipd_path)) {
    ipd_path <- "raw-data/IPDdata.txt"
  }
  if (!file.exists(ipd_path)) {
    stop("Could not locate `IPDdata.txt`.", call. = FALSE)
  }

  raw_ipd <- utils::read.table(ipd_path, header = TRUE)
  progression_data <- tibble::tibble(
    t = pmax(as.numeric(raw_ipd$time), 1e-4),
    d = as.integer(raw_ipd$event)
  )

  y <- 1444L
  n <- 6380L
  p_hat <- y / n
  mu_hat <- stats::qlogis(p_hat)
  sigma_hat <- sqrt(1 / y + 1 / (n - y))

  prevalence_data <- tibble::tibble(
    y = y,
    n = n,
    test_id = 1L,
    mu_hat = mu_hat,
    sigma_hat = sigma_hat
  )

  list(
    progression_data = progression_data,
    prevalence_data = prevalence_data
  )
}

#' Load cleaned individual-level PREDICT-TB migrant cohort data
#'
#' Loads the individual-level migrant sub-cohort from `cleaned_migrant_predict_data.RData`,
#' containing baseline diagnostic test positivity (`x`), follow-up time (`t`),
#' active TB progression indicator (`d`), centered age (`age_z`), and
#' `ethnicity`.
#'
#' @return A list with `progression_data` (tibble with individual `t`, `d`, `x`,
#'   `age`, `age_z`, `ethnicity`) and `prevalence_data` (aggregate summary of
#'   test positivity in the cohort).
#' @export
load_predict_migrant <- function() {
  rdata_path <- system.file(
    "extdata",
    "cleaned_migrant_predict_data.RData",
    package = "EvidSynthTB"
  )
  if (!nzchar(rdata_path) || !file.exists(rdata_path)) {
    rdata_path <- "data input/cleaned_migrant_predict_data.RData"
  }
  if (!file.exists(rdata_path)) {
    stop("Could not locate `cleaned_migrant_predict_data.RData`.", call. = FALSE)
  }

  env <- new.env(parent = emptyenv())
  load(rdata_path, envir = env)
  dat_m <- env$dat_m
  dat_m <- dat_m[!is.na(dat_m$age) & dat_m$ethnicity != "", , drop = FALSE]

  age_num <- as.numeric(dat_m$age)
  progression_data <- tibble::tibble(
    t = pmax(as.numeric(dat_m$time), 1e-4),
    d = as.integer(dat_m$status),
    x = as.integer(dat_m$pos),
    x1 = as.integer(dat_m$pos),
    age = age_num,
    age_z = as.numeric(scale(age_num)),
    ethnicity = droplevels(as.factor(dat_m$ethnicity))
  )

  y <- sum(progression_data$x, na.rm = TRUE)
  n <- nrow(progression_data)
  p_hat <- y / n

  prevalence_data <- tibble::tibble(
    y = as.integer(y),
    n = as.integer(n),
    test_id = 1L,
    mu_hat = stats::qlogis(p_hat),
    sigma_hat = sqrt(1 / y + 1 / (n - y))
  )

  list(
    progression_data = progression_data,
    prevalence_data = prevalence_data
  )
}
