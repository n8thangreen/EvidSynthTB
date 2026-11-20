#' Locate and compile an internal CmdStan model
#'
#' @noRd
get_evidsynth_stan_model <- function(stan_file) {
  stan_path <- system.file("stan", stan_file, package = "EvidSynthTB")
  if (!nzchar(stan_path) || !file.exists(stan_path)) {
    stan_path <- file.path("inst", "stan", stan_file)
  }
  if (!file.exists(stan_path)) {
    stop(sprintf("Could not find Stan file `%s`.", stan_file), call. = FALSE)
  }

  inc_path <- dirname(stan_path)
  inc_sub <- file.path(inc_path, "include")
  cache_dir <- tools::R_user_dir("EvidSynthTB", which = "cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  cmdstanr::cmdstan_model(
    stan_file = stan_path,
    include_paths = c(inc_path, inc_sub),
    dir = cache_dir,
    quiet = TRUE
  )
}

#' Build a numeric covariate design matrix (without intercept)
#'
#' @noRd
build_design_matrix <- function(df, covariates) {
  if (is.null(covariates) || length(covariates) == 0L || nrow(df) == 0L) {
    return(matrix(0.0, nrow = nrow(df), ncol = 0L))
  }
  missing_cols <- setdiff(covariates, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("Covariate column(s) not found: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  sub_df <- as.data.frame(df[, covariates, drop = FALSE])
  form <- stats::as.formula(paste("~", paste(covariates, collapse = " + ")))
  mm <- stats::model.matrix(form, data = sub_df)
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  storage.mode(mm) <- "double"
  mm
}

#' Fit Bayesian Multi-Parameter Evidence Synthesis Model for LTBI and TB Progression
#'
#' Links latent tuberculosis infection (LTBI) diagnostic test / survey data with
#' active TB progression data using a parametric mixture cure model in CmdStan.
#' Supports both **published aggregate / surveillance report tables** (cumulative
#' subgroup counts, UKHSA/ETS annual interval counts, and digitized Kaplan-Meier
#' life tables) and **individual-level survival / linked test data**, with
#' imperfect diagnostic tests (`sens`, `spec`) and six modular synthesis
#' strategies (`"joint"`, `"cut"`, `"independent"`, `"melding"`, `"elastic"`,
#' `"tempered"`).
#'
#' @param prevalence_data Data frame or tibble of aggregate LTBI survey data
#'   (columns `y` or `ltbi`/`n_pos` and `n` or `pop`/`n_migrant`, plus optional
#'   `test_id` and covariates) and/or normal logit plug-in summaries (`mu_hat`
#'   and `sigma_hat`). Can be `NULL` when `progression_data` already contains
#'   group-level LTBI counts (`y` / `ltbi` / `n_pos`) or individual-level test
#'   results (`x` / `x1`, `x2`, ...).
#' @param progression_data Data frame or tibble of active TB progression data in
#'   one of four formats:
#'   \enumerate{
#'     \item **Published cumulative subgroup summary table:** columns `n` (or
#'       `pop` / `n_migrant`), `d` (or `tb` / `n_tb`), optional `y` (or `ltbi` /
#'       `n_pos`), optional `t_end` (follow-up duration, defaulting to
#'       `t_horizon`), and optional subgroup covariates.
#'     \item **UKHSA / ETS annual or interval surveillance table:** columns
#'       `t_start` (`t0`), `t_end` (`t1`), `d` (`tb` / `n_tb` / `n.event`), and
#'       `n` (`pop` / `n_migrant`, initial cohort denominator).
#'     \item **Published Kaplan-Meier life table:** columns `t_start`/`t_end`
#'       (or `time`), `n_risk` (`n.risk`), and `d` (`n.event`).
#'     \item **Individual-level survival data (IPD):** columns `t` (follow-up
#'       time `> 0`) and `d` (`1` = active TB progression, `0` = right-censored),
#'       plus optional individual diagnostic test columns (`x` or `x1`, `x2`, ...).
#'   }
#' @param model Character. Synthesis model strategy:
#'   \describe{
#'     \item{`"joint"`}{Full joint synthesis sharing baseline LTBI prevalence
#'       (`alpha_surv = alpha_diag`).}
#'     \item{`"cut"`}{True two-stage Bayesian cut model (Plummer 2015): fits
#'       Stage 1 on the LTBI diagnostic/survey data alone, then propagates
#'       posterior draws of prevalence and test accuracy into Stage 2 survival
#'       inference with zero feedback from survival to prevalence.}
#'     \item{`"independent"`}{Disjoint modules with separate independent
#'       prevalence parameters `prev_diag` and `prev_cf`.}
#'     \item{`"melding"`}{Two-stage normal portal approximation (`mu_hat`,
#'       `sigma_hat` on `logit(pi)`) combined with the progression survival
#'       likelihood.}
#'     \item{`"elastic"`}{Hierarchical coupling where
#'       `logit(prev_cf) ~ Normal(logit(prev_diag), tau_couple)`.}
#'     \item{`"tempered"`}{Power-likelihood synthesis weighting the progression
#'       survival likelihood by `w_temper` in `(0, 1]`.}
#'   }
#' @param coupling Character. How diagnostic test results and progression times
#'   are coupled when individual-level linked test data are provided:
#'   `"population"` (modular population-level prevalence link) or
#'   `"latent_class"` (individual-level marginalization over each person's latent
#'   LTBI state `Z_i`).
#' @param cohort_type Character. `"all"` if all screened individuals are
#'   followed for progression (uncured fraction is `prev_cf`), or `"pos_only"`
#'   if only test-positive individuals are followed (uncured fraction is the
#'   positive predictive value, `ppv_cf`).
#' @param surv_dist Character. Parametric uncured survival family:
#'   `"gompertz"`, `"weibull"`, or `"exponential"`.
#' @param imperfect_test Logical. If `TRUE`, estimates test sensitivity (`sens`)
#'   and specificity (`spec`) with `Beta(a_sens, b_sens)` and
#'   `Beta(a_spec, b_spec)` priors. If `FALSE`, fixes them to `sens_fixed` and
#'   `spec_fixed`.
#' @param covariates Optional character vector of covariate column names to
#'   include in the linear predictor for `logit(pi)`.
#' @param t_horizon Numeric (`> 0`). Default follow-up duration (in years) for
#'   published cumulative subgroup summary tables when `t_end` is not explicitly
#'   provided as a column. Default is `5`.
#' @param priors Named list of prior hyperparameters overriding defaults (`mu_alpha`,
#'   `sigma_alpha`, `sigma_beta`, `mu_shape`, `sigma_shape`, `a_lambda`,
#'   `b_lambda`, `a_sens`, `b_sens`, `a_spec`, `b_spec`, `sens_fixed`,
#'   `spec_fixed`).
#' @param tau_couple Numeric (`> 0`). Coupling standard deviation on the logit
#'   scale when `model = "elastic"`. Default is `0.5`.
#' @param w_temper Numeric in `[0, 1]`. Power-likelihood weight on the survival
#'   module when `model = "tempered"`. Default is `0.5`.
#' @param t_lim Integer (`>= 1`). Maximum integer time horizon for posterior
#'   predicted survival curves `S_pred[1:t_lim]`. Default is `20`.
#' @param chains Integer. Number of MCMC chains. Default is `2`.
#' @param n_iter Integer. Total iterations per chain (or post-warmup draws +
#'   `n_burnin`). Default is `1000`.
#' @param n_burnin Integer. Warmup iterations per chain. Default is `500`.
#' @param n_thin Integer. Thinning interval. Default is `1`.
#' @param n_cut_subsamples Integer. Number of Stage 1 posterior draws to
#'   propagate in Stage 2 conditional HMC when `model = "cut"`. Default is `25`.
#' @param seed Optional integer random seed.
#' @param refresh Integer. CmdStan output refresh rate (`0` suppresses progress
#'   printing). Default is `0`.
#'
#' @return An S3 list object of class `"evidsynth_fit"` containing:
#' \describe{
#'   \item{`draws`}{A `posterior::draws_df` data frame of posterior samples for
#'     `prev_cf`, `prev_diag`, `ppv_cf`, `sens`, `spec`, `lambda`, `shape`,
#'     `S_pred[1:t_lim]`, `S_pop_pred[1:t_lim]`, and `log_lik_surv`.}
#'   \item{`summary`}{Posterior summary table produced by
#'     `posterior::summarise_draws()`.}
#'   \item{`fit`}{Underlying `CmdStanMCMC` object (or Stage 1/2 fit objects for
#'     `"cut"`).}
#'   \item{`data`}{Constructed Stan data list.}
#'   \item{`model`}{Character name of the synthesis strategy.}
#'   \item{`surv_dist`}{Character name of the uncured survival distribution.}
#' }
#' @export
evidsynth_fit <- function(
  prevalence_data = NULL,
  progression_data,
  model = c("cut", "independent", "joint", "melding", "elastic", "tempered"),
  coupling = c("population", "latent_class"),
  cohort_type = c("all", "pos_only"),
  surv_dist = c("gompertz", "weibull", "exponential"),
  imperfect_test = FALSE,
  covariates = NULL,
  t_horizon = 5,
  priors = list(),
  tau_couple = 0.5,
  w_temper = 0.5,
  t_lim = 20,
  chains = 2,
  n_iter = 1000,
  n_burnin = 500,
  n_thin = 1,
  n_cut_subsamples = 25,
  seed = 123,
  refresh = 0
) {
  model <- match.arg(model)
  coupling <- match.arg(coupling)
  cohort_type <- match.arg(cohort_type)
  surv_dist <- match.arg(surv_dist)

  model_type_map <- c(
    joint = 1L,
    independent = 2L,
    cut = 3L,
    melding = 4L,
    elastic = 5L,
    tempered = 6L
  )
  surv_dist_map <- c(
    gompertz = 1L,
    weibull = 2L,
    exponential = 3L
  )

  pnames <- names(progression_data)
  is_agg_surv <- any(
    c("n", "pop", "n_migrant", "n_screened", "n_risk", "n.risk", "t_start", "t0") %in% pnames
  )

  # 1. If prevalence_data is NULL and progression_data is an aggregate table
  #    with LTBI test counts (y / ltbi / n_pos), extract prevalence_data
  if (is.null(prevalence_data) && is_agg_surv) {
    y_col <- intersect(c("y", "ltbi", "n_pos"), pnames)
    n_col <- intersect(c("n", "pop", "n_migrant", "n_screened"), pnames)
    if (length(y_col) > 0L && length(n_col) > 0L) {
      prevalence_data <- progression_data
      prevalence_data$y <- as.integer(progression_data[[y_col[1]]])
      prevalence_data$n <- as.integer(progression_data[[n_col[1]]])
    }
  }

  # 2. Parse progression_data: Individual IPD vs Published Aggregate / Interval / KM Table
  if (!is_agg_surv) {
    # Individual-level survival data (t, d)
    N_surv <- nrow(progression_data)
    t_vec <- pmax(as.numeric(progression_data$t), 1e-5)
    d_vec <- as.numeric(progression_data$d)
    X_surv <- build_design_matrix(progression_data, covariates)
    P_cov <- ncol(X_surv)

    test_cols <- grep("^x[0-9]+$", pnames, value = TRUE)
    if (length(test_cols) == 0L && "x" %in% pnames) {
      test_cols <- "x"
    }
    if (length(test_cols) > 0L) {
      K_tests <- length(test_cols)
      is_linked <- 1L
      x_linked <- as.matrix(progression_data[, test_cols, drop = FALSE])
      storage.mode(x_linked) <- "integer"
    } else {
      K_tests <- 1L
      is_linked <- 0L
      x_linked <- matrix(0L, nrow = N_surv, ncol = 1L)
    }

    G_surv <- 0L
    t0_agg <- numeric(0)
    t1_agg <- numeric(0)
    d_agg <- integer(0)
    n_risk_agg <- integer(0)
    is_conditional_risk <- integer(0)
    X_surv_agg <- matrix(0.0, nrow = 0L, ncol = P_cov)
  } else {
    # Published Aggregate / Interval / KM Life Table
    N_surv <- 0L
    t_vec <- numeric(0)
    d_vec <- numeric(0)
    K_tests <- 1L
    is_linked <- 0L
    x_linked <- matrix(0L, nrow = 0L, ncol = 1L)

    prog_df <- as.data.frame(progression_data)

    # Check if raw KM life table with a single 'time' column and 'n.risk'/'n_risk'
    if (
      "time" %in% names(prog_df) &&
      !any(c("t_start", "t0", "t_end", "t1") %in% names(prog_df)) &&
      any(c("n.risk", "n_risk") %in% names(prog_df))
    ) {
      nr <- nrow(prog_df)
      if (nr < 2L) {
        stop("Kaplan-Meier table with `time` must have at least 2 rows.", call. = FALSE)
      }
      t0_raw <- as.numeric(prog_df$time)
      t1_raw <- c(t0_raw[-1], t0_raw[nr] + median(diff(t0_raw)))
      prog_df$t_start <- pmax(t0_raw, 0)
      prog_df$t_end <- pmax(t1_raw, prog_df$t_start + 1e-4)
    }

    G_surv <- nrow(prog_df)
    X_surv_agg <- build_design_matrix(prog_df, covariates)
    P_cov <- ncol(X_surv_agg)
    X_surv <- matrix(0.0, nrow = 0L, ncol = P_cov)

    # Interval start t0_agg and end t1_agg
    t0_col <- intersect(c("t_start", "t0"), names(prog_df))
    t1_col <- intersect(c("t_end", "t1", "time"), names(prog_df))
    t0_agg <- if (length(t0_col) > 0L) pmax(as.numeric(prog_df[[t0_col[1]]]), 0.0) else rep(0.0, G_surv)
    t1_agg <- if (length(t1_col) > 0L) pmax(as.numeric(prog_df[[t1_col[1]]]), t0_agg + 1e-4) else rep(as.numeric(t_horizon), G_surv)

    # Event counts d_agg
    d_col <- intersect(c("d", "tb", "n_tb", "n.event", "n_event", "events"), names(prog_df))
    if (length(d_col) == 0L) {
      stop("Aggregate `progression_data` must include an active TB event count column (`d`, `tb`, `n_tb`, or `n.event`).", call. = FALSE)
    }
    d_agg <- as.integer(round(prog_df[[d_col[1]]]))

    # Denominator n_risk_agg and conditional vs initial cohort flag
    risk_col <- intersect(c("n_risk", "n.risk"), names(prog_df))
    cohort_col <- intersect(c("n", "pop", "n_migrant", "n_screened"), names(prog_df))
    if (length(risk_col) > 0L) {
      n_risk_agg <- pmax(as.integer(round(prog_df[[risk_col[1]]])), d_agg, 1L)
      is_conditional_risk <- rep(1L, G_surv)
    } else if (length(cohort_col) > 0L) {
      n_risk_agg <- pmax(as.integer(round(prog_df[[cohort_col[1]]])), d_agg, 1L)
      is_conditional_risk <- rep(0L, G_surv)
    } else {
      stop("Aggregate `progression_data` must include a denominator column (`n`, `pop`, `n_migrant`, or `n_risk`/`n.risk`).", call. = FALSE)
    }
  }

  # 3. Parse prevalence_data (aggregate binomial y/n and/or normal mu_hat/sigma_hat)
  M_agg <- 0L
  y_agg <- integer(0)
  n_agg <- integer(0)
  test_id_agg <- integer(0)
  X_agg <- matrix(0.0, nrow = 0L, ncol = P_cov)
  M_plugin <- 0L
  mu_hat <- numeric(0)
  sigma_hat <- numeric(0)

  if (!is.null(prevalence_data) && nrow(prevalence_data) > 0L) {
    prev_names <- names(prevalence_data)
    y_col_p <- intersect(c("y", "ltbi", "n_pos"), prev_names)
    n_col_p <- intersect(c("n", "pop", "n_migrant", "n_screened"), prev_names)

    if (length(y_col_p) > 0L && length(n_col_p) > 0L && model != "melding") {
      M_agg <- nrow(prevalence_data)
      y_agg <- as.integer(round(prevalence_data[[y_col_p[1]]]))
      n_agg <- pmax(as.integer(round(prevalence_data[[n_col_p[1]]])), y_agg, 1L)
      if ("test_id" %in% prev_names) {
        test_id_agg <- as.integer(prevalence_data$test_id)
      } else {
        test_id_agg <- rep(1L, M_agg)
      }
      if (P_cov > 0L && all(covariates %in% prev_names)) {
        X_agg <- build_design_matrix(prevalence_data, covariates)
      } else {
        X_agg <- matrix(0.0, nrow = M_agg, ncol = P_cov)
      }
    } else if (all(c("mu_hat", "sigma_hat") %in% prev_names)) {
      M_plugin <- nrow(prevalence_data)
      mu_hat <- as.numeric(prevalence_data$mu_hat)
      sigma_hat <- as.numeric(prevalence_data$sigma_hat)
    } else if (length(y_col_p) > 0L && length(n_col_p) > 0L && model == "melding") {
      # Automatically compute empirical logit normal summary if melding is requested on y/n table
      y_tmp <- as.numeric(prevalence_data[[y_col_p[1]]])
      n_tmp <- as.numeric(prevalence_data[[n_col_p[1]]])
      p_tmp <- pmax(pmin((y_tmp + 0.5) / (n_tmp + 1.0), 0.99), 0.01)
      M_plugin <- nrow(prevalence_data)
      mu_hat <- stats::qlogis(p_tmp)
      sigma_hat <- sqrt(1 / (y_tmp + 0.5) + 1 / (n_tmp - y_tmp + 0.5))
    }
  }

  # 4. Merge default priors with user-supplied priors
  default_priors <- list(
    mu_alpha = stats::qlogis(0.30),
    sigma_alpha = 1.5,
    sigma_beta = 1.5,
    mu_shape = if (surv_dist == "weibull") 0.0 else 0.1,
    sigma_shape = 0.5,
    a_lambda = 1.0,
    b_lambda = 10.0,
    a_sens = rep(85.0, K_tests),
    b_sens = rep(15.0, K_tests),
    a_spec = rep(95.0, K_tests),
    b_spec = rep(5.0, K_tests),
    sens_fixed = rep(1.0, K_tests),
    spec_fixed = rep(1.0, K_tests)
  )
  p <- utils::modifyList(default_priors, priors)
  if (length(p$a_sens) == 1L && K_tests > 1L) p$a_sens <- rep(p$a_sens, K_tests)
  if (length(p$b_sens) == 1L && K_tests > 1L) p$b_sens <- rep(p$b_sens, K_tests)
  if (length(p$a_spec) == 1L && K_tests > 1L) p$a_spec <- rep(p$a_spec, K_tests)
  if (length(p$b_spec) == 1L && K_tests > 1L) p$b_spec <- rep(p$b_spec, K_tests)
  if (length(p$sens_fixed) == 1L && K_tests > 1L) p$sens_fixed <- rep(p$sens_fixed, K_tests)
  if (length(p$spec_fixed) == 1L && K_tests > 1L) p$spec_fixed <- rep(p$spec_fixed, K_tests)

  n_sampling <- max(as.integer(n_iter - n_burnin), 200L)
  n_warmup <- max(as.integer(n_burnin), 100L)

  dat_input <- list(
    model_type = unname(model_type_map[[model]]),
    coupling = if (coupling == "latent_class") 2L else 1L,
    cohort_type = if (cohort_type == "pos_only") 2L else 1L,
    surv_dist = unname(surv_dist_map[[surv_dist]]),
    imperfect_test = as.integer(imperfect_test),
    t_lim = as.integer(t_lim),
    P_cov = P_cov,
    N_surv = N_surv,
    t = array(t_vec, dim = N_surv),
    d = array(d_vec, dim = N_surv),
    X_surv = X_surv,
    G_surv = G_surv,
    t0_agg = array(t0_agg, dim = G_surv),
    t1_agg = array(t1_agg, dim = G_surv),
    d_agg = array(d_agg, dim = G_surv),
    n_risk_agg = array(n_risk_agg, dim = G_surv),
    is_conditional_risk = array(is_conditional_risk, dim = G_surv),
    X_surv_agg = X_surv_agg,
    K_tests = K_tests,
    is_linked = is_linked,
    x_linked = x_linked,
    M_agg = M_agg,
    y_agg = array(y_agg, dim = M_agg),
    n_agg = array(n_agg, dim = M_agg),
    test_id_agg = array(test_id_agg, dim = M_agg),
    X_agg = X_agg,
    N_ind = 0L,
    X_ind = matrix(0.0, nrow = 0L, ncol = P_cov),
    x_ind = array(0L, dim = c(0L, K_tests)),
    M_plugin = M_plugin,
    mu_hat = array(mu_hat, dim = M_plugin),
    sigma_hat = array(sigma_hat, dim = M_plugin),
    fixed_alpha_diag = p$mu_alpha,
    fixed_beta_cov = array(rep(0.0, P_cov), dim = P_cov),
    fixed_sens = array(p$sens_fixed, dim = K_tests),
    fixed_spec = array(p$spec_fixed, dim = K_tests),
    sens_fixed = array(p$sens_fixed, dim = K_tests),
    spec_fixed = array(p$spec_fixed, dim = K_tests),
    a_sens = array(p$a_sens, dim = K_tests),
    b_sens = array(p$b_sens, dim = K_tests),
    a_spec = array(p$a_spec, dim = K_tests),
    b_spec = array(p$b_spec, dim = K_tests),
    mu_alpha = p$mu_alpha,
    sigma_alpha = p$sigma_alpha,
    sigma_beta = p$sigma_beta,
    mu_shape = p$mu_shape,
    sigma_shape = p$sigma_shape,
    a_lambda = p$a_lambda,
    b_lambda = p$b_lambda,
    tau_couple = tau_couple,
    w_temper = w_temper
  )

  mod_modular <- get_evidsynth_stan_model("evidsynth_modular.stan")

  # --- Case A: True Two-Stage Bayesian Cut Model ---
  if (model == "cut") {
    if (M_agg > 0L || is_linked == 1L) {
      mod_ltbi <- get_evidsynth_stan_model("ltbi_module.stan")
      ltbi_dat <- list(
        K_tests = K_tests,
        P_cov = P_cov,
        M_agg = M_agg,
        y_agg = array(y_agg, dim = M_agg),
        n_agg = array(n_agg, dim = M_agg),
        test_id_agg = array(test_id_agg, dim = M_agg),
        X_agg = X_agg,
        N_ind = if (is_linked == 1L) N_surv else 0L,
        X_ind = if (is_linked == 1L) X_surv else matrix(0.0, nrow = 0L, ncol = P_cov),
        x_ind = if (is_linked == 1L) x_linked else array(0L, dim = c(0L, K_tests)),
        imperfect_test = as.integer(imperfect_test),
        sens_fixed = array(p$sens_fixed, dim = K_tests),
        spec_fixed = array(p$spec_fixed, dim = K_tests),
        a_sens = array(p$a_sens, dim = K_tests),
        b_sens = array(p$b_sens, dim = K_tests),
        a_spec = array(p$a_spec, dim = K_tests),
        b_spec = array(p$b_spec, dim = K_tests),
        mu_alpha = p$mu_alpha,
        sigma_alpha = p$sigma_alpha,
        sigma_beta = p$sigma_beta
      )
      fit_stage1 <- mod_ltbi$sample(
        data = ltbi_dat,
        chains = chains,
        iter_warmup = n_warmup,
        iter_sampling = n_sampling,
        thin = n_thin,
        seed = seed,
        refresh = refresh,
        show_messages = FALSE,
        show_exceptions = FALSE
      )
      draws_stage1 <- posterior::as_draws_df(fit_stage1$draws())
    } else if (M_plugin > 0L) {
      prec_prior <- 1 / (p$sigma_alpha^2)
      prec_data <- sum(1 / (sigma_hat^2))
      post_var <- 1 / (prec_prior + prec_data)
      post_mean <- post_var * (p$mu_alpha * prec_prior + sum(mu_hat / (sigma_hat^2)))
      n_total_draws <- chains * floor(n_sampling / n_thin)
      if (!is.null(seed)) set.seed(seed)
      alpha_draws <- stats::rnorm(n_total_draws, mean = post_mean, sd = sqrt(post_var))
      draws_stage1 <- tibble::tibble(
        alpha_diag = alpha_draws,
        prev_diag = stats::plogis(alpha_draws)
      )
      if (P_cov > 0L) {
        for (j in seq_len(P_cov)) {
          draws_stage1[[sprintf("beta_cov[%d]", j)]] <- rep(0.0, n_total_draws)
        }
      }
      for (k in seq_len(K_tests)) {
        if (imperfect_test) {
          draws_stage1[[sprintf("sens[%d]", k)]] <- stats::rbeta(
            n_total_draws, p$a_sens[k], p$b_sens[k]
          )
          draws_stage1[[sprintf("spec[%d]", k)]] <- stats::rbeta(
            n_total_draws, p$a_spec[k], p$b_spec[k]
          )
        } else {
          draws_stage1[[sprintf("sens[%d]", k)]] <- rep(p$sens_fixed[k], n_total_draws)
          draws_stage1[[sprintf("spec[%d]", k)]] <- rep(p$spec_fixed[k], n_total_draws)
        }
      }
      fit_stage1 <- NULL
    } else {
      stop("`cut` model requires either `prevalence_data` or LTBI test counts in `progression_data`.", call. = FALSE)
    }

    n_s1 <- nrow(draws_stage1)
    K_sub <- min(as.integer(n_cut_subsamples), n_s1)
    sub_idx <- unique(round(seq(1, n_s1, length.out = K_sub)))
    K_sub <- length(sub_idx)
    draws_per_sub <- max(as.integer(ceiling((chains * n_sampling / n_thin) / K_sub)), 20L)

    draws_stage1_df <- as.data.frame(draws_stage1)
    combined_list <- vector("list", K_sub)
    for (m in seq_len(K_sub)) {
      idx_m <- sub_idx[m]
      dat_m <- dat_input
      dat_m$model_type <- 3L
      dat_m$fixed_alpha_diag <- as.numeric(draws_stage1_df$alpha_diag[idx_m])
      if (P_cov > 0L) {
        beta_cols <- sprintf("beta_cov[%d]", seq_len(P_cov))
        dat_m$fixed_beta_cov <- array(
          as.numeric(draws_stage1_df[idx_m, beta_cols]),
          dim = P_cov
        )
      }
      sens_cols <- sprintf("sens[%d]", seq_len(K_tests))
      spec_cols <- sprintf("spec[%d]", seq_len(K_tests))
      dat_m$fixed_sens <- array(as.numeric(draws_stage1_df[idx_m, sens_cols]), dim = K_tests)
      dat_m$fixed_spec <- array(as.numeric(draws_stage1_df[idx_m, spec_cols]), dim = K_tests)

      fit_m <- mod_modular$sample(
        data = dat_m,
        chains = 1,
        iter_warmup = min(n_warmup, 250L),
        iter_sampling = draws_per_sub,
        seed = if (!is.null(seed)) seed + m else NULL,
        refresh = 0,
        show_messages = FALSE,
        show_exceptions = FALSE
      )
      df_m <- posterior::as_draws_df(fit_m$draws())
      combined_list[[m]] <- df_m
    }

    draws_df <- dplyr::bind_rows(combined_list)
    draws_df$.chain <- 1L
    draws_df$.iteration <- seq_len(nrow(draws_df))
    draws_df$.draw <- seq_len(nrow(draws_df))
    draws_df <- posterior::as_draws_df(draws_df)
    fit_obj <- list(stage1 = fit_stage1, stage2_subsamples = K_sub)
  } else {
    # --- Case B: Single-stage CmdStan HMC (Joint, Independent, Melding, Elastic, Tempered) ---
    fit_obj <- mod_modular$sample(
      data = dat_input,
      chains = chains,
      parallel_chains = min(chains, parallel::detectCores()),
      iter_warmup = n_warmup,
      iter_sampling = n_sampling,
      thin = n_thin,
      seed = seed,
      refresh = refresh,
      show_messages = FALSE,
      show_exceptions = FALSE
    )
    draws_df <- posterior::as_draws_df(fit_obj$draws())
  }

  keep_vars <- intersect(
    c(
      "prev_diag", "prev_cf", "ppv_cf", "lambda", "shape", "sens[1]", "spec[1]",
      if (P_cov > 0L) sprintf("beta_cov[%d]", seq_len(P_cov)) else character(0)
    ),
    names(draws_df)
  )
  summary_tbl <- posterior::summarise_draws(
    posterior::subset_draws(draws_df, variable = keep_vars)
  )

  structure(
    list(
      draws = draws_df,
      summary = summary_tbl,
      fit = fit_obj,
      data = dat_input,
      model = model,
      coupling = coupling,
      cohort_type = cohort_type,
      surv_dist = surv_dist,
      imperfect_test = imperfect_test,
      stan_params = tibble::lst(n_iter, n_burnin, n_thin, chains)
    ),
    class = "evidsynth_fit"
  )
}
