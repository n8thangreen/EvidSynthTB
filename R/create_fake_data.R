#' Simulate synthetic LTBI diagnostic and TB progression survival data
#'
#' Generates synthetic data for Bayesian Multi-Parameter Evidence Synthesis
#' (MPES) linking latent tuberculosis infection (LTBI) diagnostic test data and
#' active TB progression time-to-event data. Supports imperfect diagnostic tests
#' (`sens`, `spec`), multiple diagnostic tests per person, linked or unlinked
#' cohorts, test-positive subgroup follow-up, covariates, and prior/data
#' conflict scenarios (`p_ltbi_surv != p_ltbi`).
#'
#' @param N_surv Integer. Sample size of the progression survival cohort (before
#'   any `pos_only` filtering). Default is `500`.
#' @param n_survey Integer vector. Sample sizes of aggregate LTBI prevalence
#'   surveys. Default is `c(1000)`.
#' @param p_ltbi Numeric in `(0, 1)`. True baseline LTBI prevalence in the
#'   diagnostic survey population. Default is `0.30`.
#' @param p_ltbi_surv Numeric in `(0, 1)` or `NULL`. True baseline LTBI
#'   prevalence in the progression survival cohort. Defaults to `p_ltbi` (no
#'   conflict). Set to a different value (e.g., `0.60`) to simulate selection
#'   bias / conflict between the survey and survival modules.
#' @param sens Numeric vector in `(0, 1]`. Sensitivity of each diagnostic test
#'   `k = 1, ..., K`. Default is `0.85`.
#' @param spec Numeric vector in `(0, 1]`. Specificity of each diagnostic test
#'   `k = 1, ..., K`. Default is `0.95`.
#' @param surv_dist Character. Parametric survival distribution for uncured
#'   progression times: `"gompertz"`, `"weibull"`, or `"exponential"`.
#' @param shape Numeric. Shape parameter for `"gompertz"` (can be real) or
#'   `"weibull"` (must be positive). Ignored for `"exponential"`. Default `0.1`.
#' @param rate Numeric (`> 0`). Hazard rate / scale parameter `lambda`. Default
#'   `0.1`.
#' @param t_censor Numeric (`> 0`). Maximum administrative follow-up / censoring
#'   time. Default is `10`.
#' @param cohort_type Character. Either `"all"` (all screened individuals are
#'   followed for active TB progression) or `"pos_only"` (only individuals
#'   testing positive on test 1 are followed for progression, so false positives
#'   are cured and the effective uncured fraction is the positive predictive
#'   value, PPV).
#' @param include_covariates Logical. If `TRUE`, simulates two baseline
#'   covariates (`age_z` and `high_inc_origin`) affecting `logit(p_ltbi)`.
#' @param beta_cov Numeric vector of length 2 used when
#'   `include_covariates = TRUE`. Default is `c(0.4, 0.6)`.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return A list with components:
#' \describe{
#'   \item{`progression_data`}{Tibble of individual-level survival follow-up
#'     (`t`, `d`, `z_true`, test results `x1`, ..., and optional covariates).}
#'   \item{`prevalence_data`}{Tibble of aggregate binomial survey data (`y`,
#'     `n`, `test_id`) plus normal logit plug-in summaries (`mu_hat`,
#'     `sigma_hat`).}
#'   \item{`diagnostic_ipd`}{Tibble of individual-level test results (`x1`, ...)
#'     for the progression cohort.}
#'   \item{`truth`}{List of true data-generating parameters (`p_ltbi`,
#'     `p_ltbi_surv`, `sens`, `spec`, `ppv_surv`, `surv_dist`, `shape`,
#'     `rate`, `cohort_type`).}
#' }
#' @export
create_fake_data <- function(
  N_surv = 500,
  n_survey = c(1000),
  p_ltbi = 0.30,
  p_ltbi_surv = NULL,
  sens = 0.85,
  spec = 0.95,
  surv_dist = c("gompertz", "weibull", "exponential"),
  shape = 0.1,
  rate = 0.1,
  t_censor = 10,
  cohort_type = c("all", "pos_only"),
  include_covariates = FALSE,
  beta_cov = c(0.4, 0.6),
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  surv_dist <- match.arg(surv_dist)
  cohort_type <- match.arg(cohort_type)
  if (is.null(p_ltbi_surv)) {
    p_ltbi_surv <- p_ltbi
  }

  K_tests <- length(sens)
  if (length(spec) != K_tests) {
    stop("`sens` and `spec` must have the same length.", call. = FALSE)
  }

  # 1. Covariates and individual true LTBI probabilities in survival cohort
  alpha_diag <- stats::qlogis(p_ltbi)
  alpha_surv <- stats::qlogis(p_ltbi_surv)

  if (include_covariates) {
    age_z <- stats::rnorm(N_surv, mean = 0, sd = 1)
    high_inc_origin <- stats::rbinom(N_surv, size = 1, prob = 0.5)
    X_surv <- cbind(age_z = age_z, high_inc_origin = high_inc_origin)
    eta_surv <- as.numeric(alpha_surv + X_surv %*% beta_cov)
  } else {
    X_surv <- matrix(0, nrow = N_surv, ncol = 0)
    eta_surv <- rep(alpha_surv, N_surv)
  }

  pi_surv_i <- stats::plogis(eta_surv)
  z_true <- stats::rbinom(N_surv, size = 1, prob = pi_surv_i)

  # 2. Simulate K diagnostic tests per individual in survival cohort
  x_mat <- matrix(0L, nrow = N_surv, ncol = K_tests)
  colnames(x_mat) <- paste0("x", seq_len(K_tests))
  for (k in seq_len(K_tests)) {
    prob_k <- ifelse(z_true == 1L, sens[k], 1 - spec[k])
    x_mat[, k] <- stats::rbinom(N_surv, size = 1, prob = prob_k)
  }

  # 3. Simulate uncured progression times from chosen parametric family
  if (surv_dist == "gompertz") {
    t_event <- flexsurv::rgompertz(N_surv, shape = shape, rate = rate)
  } else if (surv_dist == "weibull") {
    # Weibull PH: H(t) = rate * t^shape => U ~ Exp(1), T = (U / rate)^(1 / shape)
    u_exp <- stats::rexp(N_surv, rate = 1)
    t_event <- (u_exp / rate)^(1 / shape)
  } else {
    t_event <- stats::rexp(N_surv, rate = rate)
  }

  # Administrative / random censoring time in [0.5 * t_censor, t_censor]
  t_cens_i <- stats::runif(N_surv, min = 0.6 * t_censor, max = t_censor)

  # Cured individuals (z_true == 0) never progress to active TB
  d <- ifelse(z_true == 1L & t_event <= t_cens_i, 1L, 0L)
  t_obs <- ifelse(d == 1L, pmax(round(t_event, 4), 1e-4), round(t_cens_i, 4))

  progression_df <- tibble::tibble(
    t = t_obs,
    d = d,
    x = x_mat[, 1],
    z_true = z_true
  )
  for (k in seq_len(K_tests)) {
    progression_df[[paste0("x", k)]] <- x_mat[, k]
  }
  if (include_covariates) {
    progression_df$age_z <- X_surv[, "age_z"]
    progression_df$high_inc_origin <- X_surv[, "high_inc_origin"]
  }

  if (cohort_type == "pos_only") {
    progression_df <- progression_df[progression_df$x1 == 1L, , drop = FALSE]
  }

  # 4. Simulate aggregate LTBI survey data (Module 1)
  p_pos_survey <- p_ltbi * sens[1] + (1 - p_ltbi) * (1 - spec[1])
  y_survey <- stats::rbinom(length(n_survey), size = n_survey, prob = p_pos_survey)

  # Compute empirical logit normal approximation ( adjusting for known sens/spec
  # via Rogan-Gladen clipped to (0.005, 0.995) ) for Melding plug-in convenience
  p_hat_raw <- pmax(pmin(y_survey / n_survey, 0.99), 0.01)
  p_hat_rg <- (p_hat_raw + spec[1] - 1) / (sens[1] + spec[1] - 1)
  p_hat_rg <- pmax(pmin(p_hat_rg, 0.99), 0.01)
  se_p_raw <- sqrt(p_hat_raw * (1 - p_hat_raw) / n_survey) / (sens[1] + spec[1] - 1)
  mu_hat <- stats::qlogis(p_hat_rg)
  sigma_hat <- pmax(se_p_raw / (p_hat_rg * (1 - p_hat_rg)), 0.05)

  prevalence_df <- tibble::tibble(
    y = as.integer(y_survey),
    n = as.integer(n_survey),
    test_id = 1L,
    mu_hat = mu_hat,
    sigma_hat = sigma_hat
  )

  p_pos_surv <- p_ltbi_surv * sens[1] + (1 - p_ltbi_surv) * (1 - spec[1])
  ppv_surv <- (p_ltbi_surv * sens[1]) / p_pos_surv

  list(
    progression_data = progression_df,
    prevalence_data = prevalence_df,
    truth = list(
      p_ltbi = p_ltbi,
      p_ltbi_surv = p_ltbi_surv,
      sens = sens,
      spec = spec,
      ppv_surv = ppv_surv,
      surv_dist = surv_dist,
      shape = if (surv_dist == "exponential") 0 else shape,
      rate = rate,
      cohort_type = cohort_type,
      beta_cov = if (include_covariates) beta_cov else numeric(0)
    )
  )
}

#' Simulate published summary tables and UKHSA-style annual surveillance tables
#'
#' Generates synthetic aggregate tables matching the two most common forms of
#' published TB literature and public health surveillance reports (without
#' individual participant data):
#' \enumerate{
#'   \item `summary_table`: Subgroup-level published summary table reporting
#'     group sample size (`n`), positive LTBI test count (`y`), cumulative
#'     active TB case count (`d`) over follow-up window `[0, t_horizon]`, and
#'     subgroup covariates.
#'   \item `interval_table`: Annual surveillance report table (e.g., UKHSA / ETS
#'     reports) giving active TB case counts (`d`) in annual intervals
#'     `[t_start, t_end]` (Year 0-1, 1-2, ..., `t_horizon`) out of an initial
#'     screened / migrant cohort (`n`).
#' }
#'
#' @param n_groups Integer. Number of published demographic/clinical strata in
#'   `summary_table`. Default is `6`.
#' @param n_per_group Integer. Sample size per group in `summary_table`. Default
#'   is `500`.
#' @param n_surveillance Integer. Total denominator cohort size for the annual
#'   surveillance `interval_table`. Default is `4000`.
#' @param p_ltbi Numeric in `(0, 1)`. True baseline LTBI prevalence in the
#'   published LTBI survey table. Default is `0.30`.
#' @param p_ltbi_surv Numeric in `(0, 1)` or `NULL`. True LTBI prevalence in the
#'   surveillance progression cohort. Defaults to `p_ltbi` (set different to
#'   simulate conflict / under-reporting / case-mix shift between published
#'   surveys and routine surveillance).
#' @param sens Numeric in `(0, 1]`. Diagnostic test sensitivity. Default `0.85`.
#' @param spec Numeric in `(0, 1]`. Diagnostic test specificity. Default `0.95`.
#' @param surv_dist Character. Uncured survival family: `"gompertz"`,
#'   `"weibull"`, or `"exponential"`.
#' @param shape Numeric. Shape parameter for Gompertz or Weibull. Default `-0.3`
#'   (decaying progression hazard over years since entry).
#' @param rate Numeric (`> 0`). Initial annual progression hazard rate `lambda`.
#'   Default `0.06`.
#' @param t_horizon Integer (`>= 1`). Number of annual follow-up intervals.
#'   Default is `5`.
#' @param seed Optional integer random seed.
#'
#' @return A list with `summary_table` (cumulative subgroup counts),
#'   `interval_table` (annual interval counts by time since entry),
#'   `prevalence_data` (aggregate survey totals), and `truth` (true parameters).
#' @export
create_fake_published_data <- function(
  n_groups = 6,
  n_per_group = 500,
  n_surveillance = 4000,
  p_ltbi = 0.30,
  p_ltbi_surv = NULL,
  sens = 0.85,
  spec = 0.95,
  surv_dist = c("gompertz", "weibull", "exponential"),
  shape = -0.30,
  rate = 0.06,
  t_horizon = 5,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  surv_dist <- match.arg(surv_dist)
  if (is.null(p_ltbi_surv)) p_ltbi_surv <- p_ltbi

  surv_fn <- function(t_val) {
    if (t_val <= 0) return(1.0)
    if (surv_dist == "gompertz") {
      flexsurv::pgompertz(t_val, shape = shape, rate = rate, lower.tail = FALSE)
    } else if (surv_dist == "weibull") {
      exp(-rate * (t_val^shape))
    } else {
      exp(-rate * t_val)
    }
  }

  # 1. Cumulative subgroup summary table (e.g. published screening paper table)
  high_inc <- rep(c(0L, 1L), length.out = n_groups)
  alpha0 <- stats::qlogis(p_ltbi)
  alpha_surv0 <- stats::qlogis(p_ltbi_surv)
  beta_inc <- 0.45
  pi_diag_g <- stats::plogis(alpha0 + beta_inc * (high_inc - mean(high_inc)))
  pi_surv_g <- stats::plogis(alpha_surv0 + beta_inc * (high_inc - mean(high_inc)))

  p_pos_g <- pi_diag_g * sens + (1 - pi_diag_g) * (1 - spec)
  p_tb_cum_g <- pi_surv_g * (1 - surv_fn(t_horizon))

  n_vec <- rep(as.integer(n_per_group), n_groups)
  y_vec <- stats::rbinom(n_groups, size = n_vec, prob = p_pos_g)
  d_vec <- stats::rbinom(n_groups, size = n_vec, prob = p_tb_cum_g)

  summary_table <- tibble::tibble(
    group = seq_len(n_groups),
    high_inc_cob = high_inc,
    n = n_vec,
    y = as.integer(y_vec),
    d = as.integer(d_vec),
    t_start = 0.0,
    t_end = as.numeric(t_horizon)
  )

  # 2. Annual surveillance interval table (e.g. UKHSA / ETS cases by year since entry)
  t_starts <- seq(0, t_horizon - 1)
  t_ends <- seq(1, t_horizon)
  delta_S <- vapply(seq_along(t_starts), function(k) {
    surv_fn(t_starts[k]) - surv_fn(t_ends[k])
  }, numeric(1))

  p_interval <- p_ltbi_surv * pmax(delta_S, 1e-8)
  d_interval <- stats::rbinom(length(t_starts), size = n_surveillance, prob = p_interval)

  interval_table <- tibble::tibble(
    year_since_entry = seq_along(t_starts),
    t_start = as.numeric(t_starts),
    t_end = as.numeric(t_ends),
    n = as.integer(rep(n_surveillance, length(t_starts))),
    d = as.integer(d_interval)
  )

  # 3. Aggregate survey total table
  y_tot <- sum(y_vec)
  n_tot <- sum(n_vec)
  p_hat <- pmax(pmin(y_tot / n_tot, 0.99), 0.01)
  p_rg <- pmax(pmin((p_hat + spec - 1) / (sens + spec - 1), 0.99), 0.01)
  se_rg <- sqrt(p_hat * (1 - p_hat) / n_tot) / (sens + spec - 1)

  prevalence_data <- tibble::tibble(
    y = as.integer(y_tot),
    n = as.integer(n_tot),
    test_id = 1L,
    mu_hat = stats::qlogis(p_rg),
    sigma_hat = pmax(se_rg / (p_rg * (1 - p_rg)), 0.05)
  )

  list(
    summary_table = summary_table,
    interval_table = interval_table,
    prevalence_data = prevalence_data,
    truth = list(
      p_ltbi = p_ltbi,
      p_ltbi_surv = p_ltbi_surv,
      sens = sens,
      spec = spec,
      surv_dist = surv_dist,
      shape = if (surv_dist == "exponential") 0 else shape,
      rate = rate,
      t_horizon = t_horizon
    )
  )
}
