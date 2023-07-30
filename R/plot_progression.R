
#' Plot progression survival curve
#'
#' @param evidsynth_obj Output of model fit
#'
#' @return
#' @export
#'
plot_progression <- function(evidsynth_obj) {

  stan_output <- extract(evidsynth_obj$fit)

  plot_dat <-
    stan_output$S_pred %>%
    as.data.frame() %>%
    mutate(sim = 1:n()) %>%
    melt(id.vars = "sim",
         variable.name = "time") %>%
    mutate(time = as.numeric(gsub("V", "", time))) %>%
    group_by(time) %>%
    summarise(median = median(value),
              mean = mean(value),
              lower50 = quantile(value, probs = 0.25),
              upper50 = quantile(value, probs = 0.75),
              lower95 = quantile(value, probs = 0.025),
              upper95 = quantile(value, probs = 0.975))

  ggplot(plot_dat, aes(time, median)) +
    geom_line() +
    geom_line(aes(y = mean), linetype = 2) +
    ylab("Survival") +
    geom_ribbon(aes(x = time, ymin = lower95, ymax = upper95),
                linetype = 0,
                alpha = 0.2) +
    geom_ribbon(aes(x = time, ymin = lower50, ymax = upper50),
                linetype = 0,
                alpha = 0.2) +
    ylim(0, 1)
}
