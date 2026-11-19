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

  target_vars <- c("prev_diag", "prev_cf", "ppv_cf", "lambda", "shape", "sens[1]", "spec[1]")

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

#' Plot posterior parameter comparisons across synthesis models
#'
#' Creates a faceted density or interval plot comparing the posterior
#' distributions of key parameters (`prev_cf`, `prev_diag`, `lambda`, `shape`,
#' `sens[1]`, `spec[1]`) across synthesis strategies.
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

