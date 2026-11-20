#' Compare multiple Bayesian evidence synthesis models
#'
#' Fits multiple modular/joint synthesis strategies (e.g., `"joint"`, `"cut"`,
#' `"independent"`, `"melding"`, `"elastic"`, `"tempered"`) to the same dataset
#' and collects tidy posterior summaries and survival curves for side-by-side
#' comparison.
#'
#' @param prevalence_data Aggregate survey or plug-in prevalence data passed to
#'   [evidsynth_fit()].
#' @param progression_data Individual-level survival (and optional linked test)
#'   data passed to [evidsynth_fit()].
#' @param models Character vector of synthesis strategies to fit. Default is
#'   `c("joint", "cut", "independent", "melding", "elastic", "tempered")`.
#' @param ... Additional arguments passed to [evidsynth_fit()] (such as
#'   `imperfect_test`, `coupling`, `cohort_type`, `surv_dist`, `n_iter`,
#'   `n_burnin`, `chains`, `seed`).
#'
#' @return A list of class `"evidsynth_comparison"` with components:
#' \describe{
#'   \item{`fits`}{Named list of `"evidsynth_fit"` objects.}
#'   \item{`summary_table`}{Combined tibble of posterior summaries across all
#'     fitted models.}
#'   \item{`draws_long`}{Combined long-format tibble of key parameter posterior
#'     draws across models.}
#' }
#' @export
compare_evidsynth <- function(
  prevalence_data = NULL,
  progression_data,
  models = c("joint", "cut", "independent", "melding", "elastic", "tempered"),
  ...
) {
  fits <- stats::setNames(vector("list", length(models)), models)
  summary_list <- vector("list", length(models))
  draws_list <- vector("list", length(models))

  target_vars <- c(
    "prev_diag", "prev_cf", "ppv_cf", "lambda", "shape", "sens[1]", "spec[1]"
  )

  for (i in seq_along(models)) {
    m <- models[[i]]
    fit_m <- evidsynth_fit(
      prevalence_data = prevalence_data,
      progression_data = progression_data,
      model = m,
      ...
    )
    fits[[m]] <- fit_m

    sum_m <- fit_m$summary
    sum_m$model <- m
    summary_list[[i]] <- dplyr::relocate(sum_m, "model")

    avail_vars <- intersect(target_vars, names(fit_m$draws))
    long_m <- lapply(avail_vars, function(v) {
      tibble::tibble(
        model = m,
        variable = v,
        value = as.numeric(fit_m$draws[[v]])
      )
    }) |>
      dplyr::bind_rows()
    draws_list[[i]] <- long_m
  }

  structure(
    list(
      fits = fits,
      summary_table = dplyr::bind_rows(summary_list),
      draws_long = dplyr::bind_rows(draws_list)
    ),
    class = "evidsynth_comparison"
  )
}

#' Evaluate parameter recovery and survival curve error against ground truth
#'
#' Computes simulation recovery metrics (Bias, RMSE, Posterior SD, 95% Credible
#' Interval Coverage, and 95% CrI Width) for key parameters and Integrated
#' Absolute Error (IAE) / Integrated RMSE (IRMSE) for the uncured survival curve
#' across models in an `"evidsynth_comparison"` object.
#'
#' @param comparison_obj An `"evidsynth_comparison"` object returned by
#'   [compare_evidsynth()].
#' @param truth Either the `truth` list returned by [create_fake_data()] or a
#'   named numeric vector/list of true parameter values (e.g.,
#'   `list(prev_diag = 0.25, prev_cf = 0.25, lambda = 0.1, shape = 0.1)`).
#'
#' @return A list with two tibbles:
#' \describe{
#'   \item{`parameter_metrics`}{Tibble with columns `model`, `variable`,
#'     `truth`, `post_mean`, `post_median`, `post_sd`, `bias`, `rmse`,
#'     `cov_95`, `width_95`.}
#'   \item{`survival_metrics`}{Tibble with columns `model`, `iae_surv`, and
#'     `irmse_surv` comparing posterior median uncured survival against the true
#'     survival curve (if `rate` / `lambda` and `shape` are in `truth`).}
#' }
#' @export
evaluate_comparison <- function(comparison_obj, truth) {
  # Normalize truth names from create_fake_data()$truth if needed
  true_map <- list()
  if (!is.null(truth$p_ltbi)) true_map$prev_diag <- truth$p_ltbi
  if (!is.null(truth$p_ltbi_surv)) {
    true_map$prev_cf <- truth$p_ltbi_surv
  } else if (!is.null(truth$p_ltbi)) {
    true_map$prev_cf <- truth$p_ltbi
  }
  if (!is.null(truth$ppv_surv)) true_map$ppv_cf <- truth$ppv_surv
  if (!is.null(truth$rate)) true_map$lambda <- truth$rate
  if (!is.null(truth$lambda)) true_map$lambda <- truth$lambda
  if (!is.null(truth$shape)) true_map$shape <- truth$shape
  if (!is.null(truth$sens)) true_map[["sens[1]"]] <- truth$sens[1]
  if (!is.null(truth$spec)) true_map[["spec[1]"]] <- truth$spec[1]
  for (nm in names(truth)) {
    if (is.numeric(truth[[nm]]) && length(truth[[nm]]) == 1L && is.null(true_map[[nm]])) {
      true_map[[nm]] <- truth[[nm]]
    }
  }

  param_rows <- list()
  surv_rows <- list()

  for (m in names(comparison_obj$fits)) {
    fit_m <- comparison_obj$fits[[m]]
    draws_df <- as.data.frame(fit_m$draws)

    for (v in intersect(names(true_map), names(draws_df))) {
      val <- as.numeric(draws_df[[v]])
      tr <- as.numeric(true_map[[v]])
      q025 <- as.numeric(stats::quantile(val, probs = 0.025, na.rm = TRUE))
      q975 <- as.numeric(stats::quantile(val, probs = 0.975, na.rm = TRUE))

      param_rows[[length(param_rows) + 1L]] <- tibble::tibble(
        model = m,
        variable = v,
        truth = tr,
        post_mean = mean(val, na.rm = TRUE),
        post_median = stats::median(val, na.rm = TRUE),
        post_sd = stats::sd(val, na.rm = TRUE),
        bias = mean(val, na.rm = TRUE) - tr,
        rmse = sqrt(mean((val - tr)^2, na.rm = TRUE)),
        cov_95 = (tr >= q025 && tr <= q975),
        width_95 = q975 - q025
      )
    }

    # Survival curve error if true rate and shape are known
    s_cols <- grep("^S_pred\\[", names(draws_df), value = TRUE)
    if (length(s_cols) > 0L && !is.null(true_map$lambda)) {
      s_mat <- as.matrix(draws_df[, s_cols, drop = FALSE])
      s_hat <- apply(s_mat, 2, stats::median, na.rm = TRUE)
      times <- seq_len(ncol(s_mat))
      s_dist <- if (!is.null(truth$surv_dist)) truth$surv_dist else fit_m$surv_dist
      tr_rate <- true_map$lambda
      tr_shape <- if (!is.null(true_map$shape)) true_map$shape else 0

      if (s_dist == "gompertz") {
        s_true <- flexsurv::pgompertz(
          times, shape = tr_shape, rate = tr_rate, lower.tail = FALSE
        )
      } else if (s_dist == "weibull") {
        s_true <- exp(-tr_rate * (times^tr_shape))
      } else {
        s_true <- exp(-tr_rate * times)
      }

      surv_rows[[length(surv_rows) + 1L]] <- tibble::tibble(
        model = m,
        iae_surv = mean(abs(s_hat - s_true)),
        irmse_surv = sqrt(mean((s_hat - s_true)^2))
      )
    }
  }

  list(
    parameter_metrics = dplyr::bind_rows(param_rows),
    survival_metrics = dplyr::bind_rows(surv_rows)
  )
}

#' Diagnose inter-module conflict, feedback shift, and survival WAIC
#'
#' Quantifies information flow and potential conflict between the LTBI
#' diagnostic module and the active TB progression survival module across fitted
#' synthesis strategies, without requiring known ground-truth parameters.
#'
#' @param comparison_obj An `"evidsynth_comparison"` object returned by
#'   [compare_evidsynth()].
#'
#' @return A tibble with one row per fitted model containing:
#' \describe{
#'   \item{`model`}{Synthesis strategy name.}
#'   \item{`prev_diag_mean`, `prev_diag_sd`}{Posterior mean and SD of LTBI
#'     prevalence in the diagnostic module.}
#'   \item{`prev_cf_mean`, `prev_cf_sd`}{Posterior mean and SD of the uncured
#'     LTBI fraction in the survival module.}
#'   \item{`feedback_shift`}{Shift in `prev_diag_mean` relative to the `"cut"`
#'     model (or `"independent"` if `"cut"` was not fitted), measuring how much
#'     survival data pull the LTBI prevalence estimate.}
#'   \item{`conflict_diff_median`, `conflict_prob_pos`}{Posterior median of
#'     `prev_cf - prev_diag` and posterior probability
#'     `P(prev_cf > prev_diag | Y)` (informative in `"independent"` and
#'     `"elastic"` models).}
#'   \item{`waic_surv`, `elpd_waic`, `p_waic`}{Widely Applicable Information
#'     Criterion (WAIC) for the progression survival module computed from
#'     pointwise `log_lik_surv`.}
#' }
#' @export
diagnose_feedback <- function(comparison_obj) {
  models <- names(comparison_obj$fits)
  ref_mod <- if ("cut" %in% models) {
    "cut"
  } else if ("independent" %in% models) {
    "independent"
  } else {
    models[[1]]
  }
  ref_draws <- as.data.frame(comparison_obj$fits[[ref_mod]]$draws)
  ref_prev_mean <- mean(ref_draws$prev_diag, na.rm = TRUE)

  rows <- lapply(models, function(m) {
    draws_df <- as.data.frame(comparison_obj$fits[[m]]$draws)
    p_diag <- as.numeric(draws_df$prev_diag)
    p_cf <- as.numeric(draws_df$prev_cf)
    diff_vec <- p_cf - p_diag

    # Compute survival WAIC from log_lik_surv[1:N_surv]
    ll_cols <- grep("^log_lik_surv\\[", names(draws_df), value = TRUE)
    if (length(ll_cols) > 0L) {
      ll_mat <- as.matrix(draws_df[, ll_cols, drop = FALSE])
      # Log pointwise predictive density lppd = sum_i log(mean_s exp(ll_{s,i}))
      max_ll <- apply(ll_mat, 2, max)
      lppd_i <- max_ll + log(colMeans(exp(sweep(ll_mat, 2, max_ll, "-"))))
      p_waic_i <- apply(ll_mat, 2, stats::var)
      elpd_waic <- sum(lppd_i - p_waic_i)
      p_waic <- sum(p_waic_i)
      waic_surv <- -2 * elpd_waic
    } else {
      elpd_waic <- NA_real_
      p_waic <- NA_real_
      waic_surv <- NA_real_
    }

    tibble::tibble(
      model = m,
      prev_diag_mean = mean(p_diag, na.rm = TRUE),
      prev_diag_sd = stats::sd(p_diag, na.rm = TRUE),
      prev_cf_mean = mean(p_cf, na.rm = TRUE),
      prev_cf_sd = stats::sd(p_cf, na.rm = TRUE),
      feedback_shift = mean(p_diag, na.rm = TRUE) - ref_prev_mean,
      conflict_diff_median = stats::median(diff_vec, na.rm = TRUE),
      conflict_prob_pos = mean(diff_vec > 0, na.rm = TRUE),
      waic_surv = waic_surv,
      elpd_waic = elpd_waic,
      p_waic = p_waic
    )
  })

  dplyr::bind_rows(rows)
}

#' Plot posterior parameter comparisons across synthesis models
#'
#' Creates a faceted density plot comparing the posterior distributions of key
#' parameters (`prev_cf`, `prev_diag`, `lambda`, `shape`, `sens[1]`, `spec[1]`)
#' across synthesis strategies.
#'
#' @param comparison_obj An `"evidsynth_comparison"` object returned by
#'   [compare_evidsynth()].
#' @param variables Character vector of parameter names to display. Default is
#'   `c("prev_cf", "prev_diag", "lambda", "shape")`.
#' @param truth Optional named numeric vector or list of true parameter values
#'   (e.g. `c(prev_cf = 0.3, prev_diag = 0.3, lambda = 0.1, shape = 0.1)`) to
#'   overlay as vertical dashed reference lines.
#'
#' @return A `ggplot2` object.
#' @export
plot_model_comparison <- function(
  comparison_obj,
  variables = c("prev_cf", "prev_diag", "lambda", "shape"),
  truth = NULL
) {
  df <- comparison_obj$draws_long
  df <- df[df$variable %in% variables, , drop = FALSE]
  df$variable <- factor(df$variable, levels = variables)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$value, color = .data$model, fill = .data$model)
  ) +
    ggplot2::geom_density(alpha = 0.2) +
    ggplot2::facet_wrap(ggplot2::vars(.data$variable), scales = "free") +
    ggplot2::xlab("Parameter value") +
    ggplot2::ylab("Posterior density") +
    ggplot2::theme_bw()

  if (!is.null(truth)) {
    truth_vec <- unlist(truth)
    keep_t <- intersect(names(truth_vec), variables)
    if (length(keep_t) > 0L) {
      truth_df <- tibble::tibble(
        variable = factor(keep_t, levels = variables),
        true_val = as.numeric(truth_vec[keep_t])
      )
      p <- p +
        ggplot2::geom_vline(
          data = truth_df,
          ggplot2::aes(xintercept = .data$true_val),
          linetype = "dashed",
          linewidth = 0.7,
          inherit.aes = FALSE
        )
    }
  }

  p
}

#' Plot survival curve comparisons across synthesis models
#'
#' Plots posterior median and 95% credible intervals of the progression survival
#' curve across multiple fitted synthesis strategies from [compare_evidsynth()].
#'
#' @param comparison_obj An `"evidsynth_comparison"` object returned by
#'   [compare_evidsynth()].
#' @param type Character. `"uncured"` (`S_pred`) or `"population"`
#'   (`S_pop_pred`).
#'
#' @return A `ggplot2` object.
#' @export
plot_survival_comparison <- function(
  comparison_obj,
  type = c("uncured", "population")
) {
  type <- match.arg(type)
  prefix <- if (type == "uncured") "^S_pred\\[" else "^S_pop_pred\\["

  surv_list <- lapply(names(comparison_obj$fits), function(m) {
    draws_df <- as.data.frame(comparison_obj$fits[[m]]$draws)
    s_cols <- grep(prefix, names(draws_df), value = TRUE)
    s_mat <- as.matrix(draws_df[, s_cols, drop = FALSE])
    tibble::tibble(
      model = m,
      time = seq_len(ncol(s_mat)),
      median = apply(s_mat, 2, stats::median, na.rm = TRUE),
      lower95 = apply(s_mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
      upper95 = apply(s_mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    )
  })
  surv_df <- dplyr::bind_rows(surv_list)

  ggplot2::ggplot(
    surv_df,
    ggplot2::aes(
      x = .data$time,
      y = .data$median,
      color = .data$model,
      fill = .data$model
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower95, ymax = .data$upper95),
      alpha = 0.12,
      linetype = 0
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::xlab("Time") +
    ggplot2::ylab(
      if (type == "uncured") "Uncured survival S_u(t)" else "Population survival S(t)"
    ) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_bw()
}
