
# Fit markov melding model using Stan
# with simulated artificial ltbi and tb progression data
# without covariate (age or ethnicity)

library(rstan)
library(shinystan)
library(purrr)
library(readr)
library(dplyr)


# sample size
N <- 1000

# disease-free censoring additional time
t_offset <- 5

# prevalence of ltbi
p_ltbi <- 0.3

# progression from ltbi to active tb times
progression_dat <-
  data.frame(
    t = round(flexsurv::rgompertz(N, shape = 0.1, rate = 0.1), 3)) |>
    mutate(x = rbinom(n = N, size = 1, prob = p_ltbi),   # ltbi status
           ## observe all progression times
           d = ifelse(x == 1, 1, 0),                     # censoring status
           # d = sample(c(0,1), size = N, replace = TRUE),
           # d = ifelse(x == 0, 0, d),
           t = ifelse(x == 0, t + t_offset, t)) |>       # time
  as_tibble()

prevalence_dat <-
  data.frame(mu_hat = boot::logit(0.3),
             sigma_hat = 0.01)

# rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


dat_input <-
  list(
    N = nrow(progression_dat),
    M = nrow(prevalence_dat),
    t = as.numeric(progression_dat$t),
    d = as.numeric(progression_dat$d),
    pos = progression_dat$x,
    t_lim = 20,        # maximum time
    mu_alpha = -0.8,
    scale_alpha = 0.5,
    ## normal priors on shape
    mu_shape = 0.1,
    sigma_shape = 0.5,
    ## gamma priors on rate
    a_lambda = 0.1,
    b_lambda = 0.01,
    mu_hat = prevalence_dat$mu_hat,
    sigma_hat = prevalence_dat$sigma_hat,
    mu_cf = boot::logit(0.3),
    sigma_cf = 0.1)

params <- c(
  # "S_pred",          # generated sample
  "prev_cf",          # ltbi prevalence
  "lambda", "shape")

n_iter <- 10e3
n_burnin <- 3e1
n_thin <- 2e1  #floor((n_iter - n_burnin)/500)

###########
# run MCMC
###########

out <- stan(data = dat_input,
            pars = params,
            file = here::here("stan", "Stan_code_markov_melding.stan"),
            chains = 1,
            iter = n_iter,
            warmup = n_burnin,
            thin = n_thin,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 20))

stan_output <- extract(out)

mean(stan_output$lambda)
mean(stan_output$shape)


#######
# plot
#######

library(reshape2)
library(ggplot2)

# base R
stan_output$S_pred %>%
  as.data.frame %>%
  summarise(across(.fns = ~ mean(.x, na.rm = TRUE ))) %>%
  # summarise(across(.fns = ~ median(.x, na.rm = TRUE ))) %>%
  unlist() %>%
  plot(type = "l")


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


# LTBI prevalence

hist(stan_output$cf, breaks = 40)

plot(flexsurv::pgompertz(q = 0:20, shape = 0.5, rate = 0.01, lower.tail = FALSE),
      type = "l", col = "red")

