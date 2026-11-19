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
