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

#' Load published Kaplan-Meier interval life table from Abubakar et al. (2018)
#'
#' Loads the digitized Kaplan-Meier interval table (`KMdata.txt`, containing 69
#' time intervals with `n_risk` and `d` progression events) together with the
#' published aggregate baseline IGRA positivity count (`y = 1444` positive out
#' of `n = 6380` tested). Fitting [evidsynth_fit()] to this published life table
#' runs orders of magnitude faster than fitting to 6,380 reconstructed
#' individual rows.
#'
#' @return A list with `progression_data` (tibble with `t_start`, `t_end`,
#'   `n_risk`, `d`, `n_censored`) and `prevalence_data` (tibble with `y = 1444L`,
#'   `n = 6380L`, `test_id = 1L`, `mu_hat`, `sigma_hat`).
#' @export
load_abubakar_km <- function() {
  km_path <- system.file("extdata", "KMdata.txt", package = "EvidSynthTB")
  if (!nzchar(km_path) || !file.exists(km_path)) {
    km_path <- "raw-data/KMdata.txt"
  }
  if (!file.exists(km_path)) {
    stop("Could not locate `KMdata.txt`.", call. = FALSE)
  }

  raw_km <- utils::read.table(km_path, header = TRUE)
  nr <- nrow(raw_km)
  t0 <- pmax(as.numeric(raw_km$time), 0.0)
  t1 <- c(t0[-1], t0[nr] + 0.25)

  progression_data <- tibble::tibble(
    t_start = t0,
    t_end = pmax(t1, t0 + 1e-4),
    n_risk = as.integer(raw_km$n.risk),
    d = as.integer(raw_km$n.event),
    n_censored = as.integer(raw_km$n.censored)
  )

  y <- 1444L
  n <- 6380L
  p_hat <- y / n
  prevalence_data <- tibble::tibble(
    y = y,
    n = n,
    test_id = 1L,
    mu_hat = stats::qlogis(p_hat),
    sigma_hat = sqrt(1 / y + 1 / (n - y))
  )

  list(
    progression_data = progression_data,
    prevalence_data = prevalence_data
  )
}

#' Load published-style aggregated PREDICT-TB subgroup contingency table
#'
#' Loads the grouped contingency table (`aggregated_data.csv`) summarizing
#' sample size (`pop`), positive LTBI test counts (`ltbi`), and cumulative
#' active TB progression counts (`tb`) stratified by `age_grp`, `sex`,
#' `ethnicity`, TB incidence in country of birth (`inc_cob_participant2`),
#' years since UK entry (`yearssinceentry_grp`), prior BCG vaccination
#' (`prevbcg`), and `reasonforscreening` (`"Contact"` vs `"Migrant"`).
#'
#' @param collapse_by Optional character vector of grouping columns to collapse
#'   over (summing `pop`, `ltbi`, and `tb`). Defaults to
#'   `c("reasonforscreening", "inc_cob_participant2", "age_grp")`. Set to `NULL`
#'   to return all raw cells.
#' @param t_horizon Numeric (`> 0`). Cumulative follow-up horizon in years
#'   assigned to `t_end`. Default is `5`.
#'
#' @return A list with `summary_table` (tibble with `n`, `y`, `d`, `t_start`,
#'   `t_end`, and grouping covariates) and `prevalence_data` (same table ready
#'   for [evidsynth_fit()]).
#' @export
load_predict_aggregated <- function(
  collapse_by = c("reasonforscreening", "inc_cob_participant2", "age_grp"),
  t_horizon = 5
) {
  agg_path <- system.file("extdata", "aggregated_data.csv", package = "EvidSynthTB")
  if (!nzchar(agg_path) || !file.exists(agg_path)) {
    agg_path <- "raw-data/aggregated_data.csv"
  }
  if (!file.exists(agg_path)) {
    stop("Could not locate `aggregated_data.csv`.", call. = FALSE)
  }

  raw_df <- utils::read.csv(agg_path, stringsAsFactors = FALSE)
  if ("X" %in% names(raw_df)) {
    raw_df$X <- NULL
  }

  for (col in intersect(
    c("age_grp", "sex", "ethnicity", "inc_cob_participant2", "yearssinceentry_grp", "prevbcg", "reasonforscreening"),
    names(raw_df)
  )) {
    val <- as.character(raw_df[[col]])
    val[is.na(val) | val == "" | val == "NA"] <- "Missing"
    raw_df[[col]] <- as.factor(val)
  }

  if (!is.null(collapse_by) && length(collapse_by) > 0L) {
    summary_df <- raw_df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(collapse_by))) |>
      dplyr::summarise(
        n = as.integer(sum(.data$pop, na.rm = TRUE)),
        y = as.integer(sum(.data$ltbi, na.rm = TRUE)),
        d = as.integer(sum(.data$tb, na.rm = TRUE)),
        .groups = "drop"
      )
  } else {
    summary_df <- tibble::as_tibble(raw_df) |>
      dplyr::mutate(
        n = as.integer(.data$pop),
        y = as.integer(.data$ltbi),
        d = as.integer(.data$tb)
      )
  }

  summary_df <- summary_df |>
    dplyr::filter(.data$n > 0L) |>
    dplyr::mutate(
      t_start = 0.0,
      t_end = as.numeric(t_horizon)
    )

  list(
    summary_table = summary_df,
    prevalence_data = summary_df
  )
}
