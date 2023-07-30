
#' Evidence synthesis model fit
#'
#' @param prevalence_data
#' @param progression_data
#' @param cut Cut-point model to stop feedback; logical
#' @param n_iter Number of iterations
#' @param n_burnin Number of iterations in burn-in
#' @param n_thin Number of iterations to thin
#'
#' @return List
#' @export
#'
evidsynth_fit <- function(prevalence_data,
                          progression_data,
                          cut = TRUE,
                          n_iter = 2e3,
                          n_burnin = 100,
                          n_thin = 10) {

  # rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  dat_input <-
    list(
      N = nrow(progression_data),
      M = nrow(prevalence_data),
      t = as.numeric(progression_data$t),
      d = as.numeric(progression_data$d),
      pos = progression_data$x,
      t_lim = 20,              # maximum time
      ## normal prior on shape
      mu_shape = 0.1,
      sigma_shape = 0.05,
      ## gamma prior on rate
      a_lambda = 0.2,
      b_lambda = 0.2,
      mu_hat = prevalence_data$mu_hat,
      sigma_hat = prevalence_data$sigma_hat,
      # mu_hat = array(prevalence_data$mu_hat, 1),      # for single data point
      # sigma_hat = array(prevalence_data$sigma_hat, 1),
      ## normal prior on linear transformed cure fraction
      mu_cf = boot::logit(0.3),
      sigma_cf = 0.2)

  params <- c(
    "S_pred",           # generated values
    "prev_cf",          # ltbi prevalence
    "lambda", "shape")

  ###########
  # run MCMC
  ###########

  if (cut) model_nm <- "Stan_code_markov_melding.stan"

  out <- stan(data = dat_input,
              pars = params,
              model_name = model_nm,
              # init = ,
              file = here::here("stan", model_nm),
              chains = 1,
              iter = n_iter,
              warmup = n_burnin,
              thin = n_thin,
              control = list(adapt_delta = 0.99,
                             max_treedepth = 20))

  list(fit = out,
       data = dat_input,
       params = params,
       stan_params =
         tibble::lst(n_iter,
                     n_burnin,
                     n_thin))
}
