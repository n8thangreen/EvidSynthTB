
# simulates artificial ltbi and tb data
# without covariate age or ethnicity
# and with imperfect diagnostic test

library(rstan)
library(shinystan)
library(purrr)
library(readr)
library(dplyr)

# rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

N <- 1000

rdat <-
  data.frame(
    t = round(flexsurv::rgompertz(N, shape = 0.5, rate = 0.1), 3)) |>
    mutate(x = rbinom(n = N, size = 1, prob = 0.3),
           ## observe all progression times
           ## no false negatives
           d = ifelse(x == 1, 1, 0),
           # d = sample(c(0,1), size = N, replace = TRUE),
           # d = ifelse(x == 0, 0, d),
           ## censored times after event times
           t = ifelse(x == 0, t + 5, t)) |>
  as_tibble()

rdat


dat_input <-
  list(
    N = nrow(rdat),
    N_obs = sum(rdat$d),
    N_cens = sum(1 - rdat$d),
    t_lim = 20,
    mu_shape = 0.5,
    sigma_shape = 1,
    a_lambda = 3,
    b_lambda = 30,
    obs_idx = which(rdat$d == 1),
    cens_idx = which(rdat$d == 0),
    t = as.numeric(rdat$t),
    d = as.numeric(rdat$d),
    x = rdat$x)

params <- c(
  "S_pred",
  "sens", "spec",
  "prev_cf", "prev_diag", "prev_mean",
  "lambda", "shape")

n_iter <- 10e3
n_burnin <- 3e1
n_thin <- 2e1 #floor((n_iter - n_burnin)/500)


###########
# run MCMC
###########

out <- stan(data = dat_input,
            pars = params,
            file = here::here("BUGS", "Stan_code_mixture_cure_model_simple_test.stan"),
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

hist(stan_output$prev_cf, breaks = 40)
hist(stan_output$prev_diag, breaks = 40)
hist(stan_output$prev_mean, breaks = 40)
hist(1/(1 + exp(-stan_output$prev_mean)), breaks = 40)
hist(stan_output$sens, breaks = 40)
hist(stan_output$spec, breaks = 40)

plot(flexsurv::pgompertz(q = 0:20, shape = 0.5, rate = 0.1, lower.tail = FALSE),
      type = "l")

