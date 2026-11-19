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

#' Fit Bayesian Multi-Parameter Evidence Synthesis Model for LTBI and TB Progression
#'
#' Links latent tuberculosis infection (LTBI) diagnostic test / survey data with
#' active TB progression time-to-event data using a parametric mixture cure
#' model in CmdStan. Supports imperfect diagnostic tests (`sens`, `spec`),
#' multiple tests per person, linked or unlinked datasets, and six modular
#' synthesis strategies (`"joint"`, `"cut"`, `"independent"`, `"melding"`,
#' `"elastic"`, `"tempered"`).
#'
#' @param prevalence_data Data frame or tibble of aggregate LTBI survey data
#'   (columns `y` and `n`, and optional `test_id`) and/or normal logit plug-in
#'   summaries (`mu_hat` and `sigma_hat`). Can be `NULL` when individual-level
#'   diagnostic test results `x` (`x1`, `x2`, ...) are present in
#'   `progression_data`.
#' @param progression_data Data frame or tibble of time-to-event follow-up data
#'   with columns `t` (follow-up time `> 0`) and `d` (`1` = active TB
#'   progression, `0` = right-censored), plus optional individual diagnostic
#'   test columns (`x` or `x1`, `x2`, ...) and covariates.
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
#' @param covariates Optional character vector of covariate column names in
#'   `progression_data` to include in the linear predictor for `logit(pi_i)`.
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
#'     `prev_cf`, `prev_diag`, `ppv_cf`, `sens`, `spec`, `lambda`, `shape`, and
#'     `S_pred[1:t_lim]`.}
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

  # 1. Extract individual-level test columns from progression_data if present
  test_cols <- grep("^x[0-9]+$", names(progression_data), value = TRUE)
  if (length(test_cols) == 0L && "x" %in% names(progression_data)) {
    test_cols <- "x"
  }

  N_surv <- nrow(progression_data)
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

  # 2. Extract covariate design matrix X_surv
  if (!is.null(covariates) && length(covariates) > 0L) {
    X_surv <- as.matrix(progression_data[, covariates, drop = FALSE])
    storage.mode(X_surv) <- "double"
    P_cov <- ncol(X_surv)
  } else {
    P_cov <- 0L
    X_surv <- matrix(0.0, nrow = N_surv, ncol = 0L)
  }

  # 3. Parse prevalence_data (aggregate binomial y/n and/or normal mu_hat/sigma_hat)
  M_agg <- 0L
  y_agg <- integer(0)
  n_agg <- integer(0)
  test_id_agg <- integer(0)
  M_plugin <- 0L
  mu_hat <- numeric(0)
  sigma_hat <- numeric(0)

  if (!is.null(prevalence_data) && nrow(prevalence_data) > 0L) {
    if (all(c("y", "n") %in% names(prevalence_data)) && model != "melding") {
      M_agg <- nrow(prevalence_data)
      y_agg <- as.integer(prevalence_data$y)
      n_agg <- as.integer(prevalence_data$n)
      if ("test_id" %in% names(prevalence_data)) {
        test_id_agg <- as.integer(prevalence_data$test_id)
      } else {
        test_id_agg <- rep(1L, M_agg)
      }
    } else if (all(c("mu_hat", "sigma_hat") %in% names(prevalence_data))) {
      M_plugin <- nrow(prevalence_data)
      mu_hat <- as.numeric(prevalence_data$mu_hat)
      sigma_hat <- as.numeric(prevalence_data$sigma_hat)
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
    N_surv = N_surv,
    t = pmax(as.numeric(progression_data$t), 1e-5),
    d = as.numeric(progression_data$d),
    t_lim = as.integer(t_lim),
    P_cov = P_cov,
    X_surv = X_surv,
    K_tests = K_tests,
    is_linked = is_linked,
    x_linked = x_linked,
    M_agg = M_agg,
    y_agg = array(y_agg, dim = M_agg),
    n_agg = array(n_agg, dim = M_agg),
    test_id_agg = array(test_id_agg, dim = M_agg),
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
    # Stage 1: Obtain posterior draws of (alpha_diag, beta_cov, sens, spec)
    # strictly from the LTBI diagnostic/survey module (no survival feedback)
    if (M_agg > 0L || is_linked == 1L) {
      mod_ltbi <- get_evidsynth_stan_model("ltbi_module.stan")
      ltbi_dat <- list(
        K_tests = K_tests,
        M_agg = M_agg,
        y_agg = array(y_agg, dim = M_agg),
        n_agg = array(n_agg, dim = M_agg),
        test_id_agg = array(test_id_agg, dim = M_agg),
        N_ind = if (is_linked == 1L) N_surv else 0L,
        P_cov = P_cov,
        X_ind = if (is_linked == 1L) X_surv else matrix(0.0, nrow = 0L, ncol = P_cov),
        x_ind = if (is_linked == 1L) x_linked else array(0L, dim = c(0L, K_tests)),
        imperfect_test = as.integer(imperfect_test),
        sens_fixed = array(p$sens_fixed, dim = K_tests),
        spec_fixed = array(p$spec_fixed, dim = K_tests),
        a_sens = array(p$a_sens, dim = K_tests),
        b_sens = array(p$a_sens * 0 + p$b_sens, dim = K_tests),
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
      # Conjugate normal posterior for alpha_diag given normal plug-in summaries
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
      stop("`cut` model requires either `prevalence_data` or linked test data in `progression_data`.", call. = FALSE)
    }

    # Stage 2: Multiple-imputation conditional HMC over K_sub Stage 1 draws
    n_s1 <- nrow(draws_stage1)
    K_sub <- min(as.integer(n_cut_subsamples), n_s1)
    sub_idx <- unique(round(seq(1, n_s1, length.out = K_sub)))
    K_sub <- length(sub_idx)
    draws_per_sub <- max(as.integer( ceiling((chains * n_sampling / n_thin) / K_sub) ), 20L)

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
    c("prev_diag", "prev_cf", "ppv_cf", "lambda", "shape", "sens[1]", "spec[1]"),
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
