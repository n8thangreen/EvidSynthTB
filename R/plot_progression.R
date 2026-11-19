#' Plot posterior progression survival curve
#'
#' Plots the posterior median, mean, 50% credible interval, and 95% credible
#' interval of the survival curve over time from an `evidsynth_fit` object.
#'
#' @param evidsynth_obj Output of [evidsynth_fit()].
#' @param type Character. Either `"uncured"` (uncured progression survival
#'   `S_pred`, conditional on being LTBI positive) or `"population"` (overall
#'   mixture survival `S_pop_pred`, including the cured fraction `1 - prev_cf`).
#'
#' @return A `ggplot2` object.
#' @export
plot_progression <- function(
  evidsynth_obj,
  type = c("uncured", "population")
) {
  type <- match.arg(type)
  prefix <- if (type == "uncured") "^S_pred\\[" else "^S_pop_pred\\["

  draws_df <- as.data.frame(evidsynth_obj$draws)
  s_cols <- grep(prefix, names(draws_df), value = TRUE)
  if (length(s_cols) == 0L) {
    stop("No survival prediction columns found in `evidsynth_obj$draws`.", call. = FALSE)
  }

  s_mat <- as.matrix(draws_df[, s_cols, drop = FALSE])
  times <- seq_len(ncol(s_mat))

  plot_dat <- tibble::tibble(
    time = times,
    median = apply(s_mat, 2, stats::median, na.rm = TRUE),
    mean = colMeans(s_mat, na.rm = TRUE),
    lower50 = apply(s_mat, 2, stats::quantile, probs = 0.25, na.rm = TRUE),
    upper50 = apply(s_mat, 2, stats::quantile, probs = 0.75, na.rm = TRUE),
    lower95 = apply(s_mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    upper95 = apply(s_mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  )

  ggplot2::ggplot(plot_dat, ggplot2::aes(x = .data$time, y = .data$median)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower95, ymax = .data$upper95),
      alpha = 0.2
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower50, ymax = .data$upper50),
      alpha = 0.25
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_line(ggplot2::aes(y = .data$mean), linetype = 2) +
    ggplot2::ylab(
      if (type == "uncured") "Uncured survival S_u(t)" else "Population survival S(t)"
    ) +
    ggplot2::xlab("Time") +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_bw()
}
